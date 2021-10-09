#!/usr/bin/env nextflow

println("Currently using the Grandeur workflow for use with microbial Illumina MiSeq sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v0.1.20211031")
println("")

// TODO : blobtools
// TODO : mycosnp
// TODO : something for plasmids
// TODO : frp_plasmid
// TODO : pointfinder
// TODO : mubsuite
// TODO : sistr
// TODO : plasmidseeker

params.reads = workflow.launchDir + '/reads'
params.outdir = workflow.launchDir + '/grandeur'
println("The files and directory for results is " + params.outdir)

params.maxmem = Math.round(Runtime.runtime.totalMemory() / 10241024).GB
println("The maximum amount of memory used in this workflow is ${params.maxmem}")
params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

params.assemble = true
if (params.assemble) {
  Channel
    .fromFilePairs("${params.reads}/*R{1,2}*.fastq*", size: 2 )
    .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
    .ifEmpty{
      println("FATAL : No fastq or fastq.gz files were found at ${params.reads}")
      println("Set 'params.reads' to directory with paired-end reads")
      exit 1
    }
    .view { "Paired-end fastq files : ${it[0]}" }
    .into { reads_seqyclean ; reads_fastqc ; reads }
} else {
  Channel
    .empty
    .into { reads_seqyclean ; reads_fastqc ; reads }
}

params.cg_pipeline = 'true'
params.genome_sizes = workflow.projectDir + "/configs/genome_sizes.json"
Channel
  .fromPath(params.genome_sizes, type:'file')
  .view { "Genome Sizes File : $it" }
  .ifEmpty{
    println("WARNING : No file with genome sizes was found at ${params.genome_sizes}")
    println("This file is required for cg-pipeline.")
  }
  .set { genome_sizes }

params.blobtools = false // right now this isn't currently working
params.blast_db = 'blast_db'
local_blastdb = params.blobtools
              ? Channel
                  .fromPath(params.blast_db, type:'dir')
                  .ifEmpty{
                    println("No blast database was found at ${params.blast_db}")
                    println("Set 'params.blast_db' to directory with blast database")
                    exit 1
                  }
                  .view { "Local Blast Database for Blobtools : $it" }
              : Channel.empty()

params.kraken2 = false
params.kraken2_db = 'kraken2_db'
if ( params.kraken2 ) {
  Channel
    .fromPath(params.kraken2_db, type:'dir')
    .ifEmpty{
      println("No kraken2 database was found at ${params.kraken2_db}")
      println("Set 'params.kraken2_db' to directory with kraken2 database")
      exit 1
      }
    .view { "Local kraken2 database : $it" }
    .set { local_kraken2 }
  kraken2_files_empty = Channel.empty()
} else {
  kraken2_files_empty = Channel.empty()
  local_kraken2 = Channel.empty()
}

params.annotation = false
params.fastas = workflow.launchDir + '/fastas'
if ( params.annotation ) {
  Channel
    .fromPath("${params.fastas}/*{.fa,.fasta.fna}")
    .ifEmpty{
      println("FATAL : No fastas were found at ${params.fastas}")
      println("Set 'params.fastas' to directory with fastas")
      exit 1
      }
    .map { file -> tuple(file.baseName, file) }
    .view { "Fastas for serotyping and AMR gene annotation : $it" }
    .into { fastas ; fastas_mash ; fastas_quast ; fastas_prokka ; fastas_seqsero2 ; fastas_amrfinder ; fastas_serotypefinder ; fastas_kleborate ; fastas_mlst }
} else {
  Channel
    .empty()
    .into { fastas ; fastas_mash ; fastas_quast ; fastas_prokka ; fastas_seqsero2 ; fastas_amrfinder ; fastas_serotypefinder ; fastas_kleborate ; fastas_mlst }
}

params.seqyclean_minlen = 25
params.seqyclean_contaminant_file = "/Adapters_plus_PhiX_174.fasta"
params.seqyclean_options = ''
process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/seqyclean:latest'

  input:
  tuple val(sample), file(reads) from reads_seqyclean

  output:
  tuple sample, file("seqyclean/${sample}_clean_PE{1,2}.fastq.gz") into clean_reads_shovill, clean_reads_mash, clean_reads_cg, clean_reads_seqsero2, clean_reads_bwa, clean_reads_shigatyper, clean_reads_kraken2
  file("seqyclean/${sample}_clean_SE.fastq.gz")
  file("seqyclean/${sample}_clean_SummaryStatistics.tsv") into seqyclean_files, seqyclean_files_combine
  file("seqyclean/${sample}_clean_SummaryStatistics.txt")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(perc_kept) into seqyclean_perc_kept_results
  tuple sample, env(kept) into seqyclean_pairskept_results

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "seqyclean version: $(seqyclean -h | grep Version | head -n 1)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    seqyclean !{params.seqyclean_options} \
      -minlen !{params.seqyclean_minlen} \
      -qual \
      -c !{params.seqyclean_contaminant_file} \
      -1 !{reads[0]} \
      -2 !{reads[1]} \
      -o seqyclean/!{sample}_clean \
      -gz \
      2>> $err_file >> $log_file

    kept=$(cut -f 58 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
    perc_kept=$(cut -f 59 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)

    if [ -z "$kept" ] ; then kept="0" ; fi
    if [ -z "$perc_kept" ] ; then perc_kept="0" ; fi
  '''
}

seqyclean_files_combine
  .collectFile(name: "SummaryStatistics.tsv",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/seqyclean")

params.shovill = true
params.shovill_options = ''
process shovill {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  memory params.maxmem
  cpus params.maxcpus
  container 'staphb/shovill:latest'

  when:
  params.shovill

  input:
  tuple val(sample), file(reads) from clean_reads_shovill

  output:
  file("shovill/${sample}/contigs.{fa,gfa}")
  file("shovill/${sample}/shovill.{corrections,log}")
  file("shovill/${sample}/spades.fasta")
  tuple sample, file("contigs/${sample}_contigs.fa") into contigs_prokka, contigs_quast, contigs_blastn, contigs_mlst, contigs_bwa, contigs_create, contigs_kleborate, contigs_amrfinder, contigs_serotypefinder
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} contigs logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    shovill --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    memory=$(echo !{task.memory} | awk '{ print $1 }' )

    shovill !{params.shovill_options} \
      --cpu !{task.cpus} \
      --ram $memory \
      --outdir !{task.process}/!{sample} \
      --R1 !{reads[0]} \
      --R2 !{reads[1]} \
      2>> $err_file >> tee $log_file
    cp !{task.process}/!{sample}/contigs.fa contigs/!{sample}_contigs.fa
  '''
}

params.fastqc = true
params.fastqc_options = ''
process fastqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/fastqc:latest'

  when:
  params.fastqc

  input:
  tuple val(sample), file(raw) from reads_fastqc

  output:
  file("${task.process}/*.html")
  file("${task.process}/*_fastqc.zip") into fastqc_files
  tuple sample, env(raw_1) into fastqc_1_results
  tuple sample, env(raw_2) into fastqc_2_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    fastqc --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fastqc !{params.fastqc_options} \
      --outdir fastqc \
      --threads !{task.cpus} \
      !{raw} \
      2>> $err_file >> $log_file

    zipped_fastq=($(ls fastqc/*fastqc.zip) "")

    raw_1=$(unzip -p ${zipped_fastq[0]} */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=$(unzip -p fastqc/*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )

    if [ -z "$raw_1" ] ; then raw_1="0" ; fi
    if [ -z "$raw_2" ] ; then raw_2="0" ; fi
  '''
}

params.mash = true
process mash_sketch {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/mash:latest'

  when:
  params.mash

  input:
  tuple val(sample), file(reads) from clean_reads_mash

  output:
  tuple sample, file("mash/${sample}.msh") into mash_sketch_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(genome_size) into mash_genome_size_results, mash_genome_size_gc
  tuple sample, env(coverage) into mash_coverage_results

  shell:
  '''
    mkdir -p mash logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    cat !{reads} | mash sketch -m 2 -o mash/!{sample} - 2>> $err_file | tee $log_file

    genome_size=$(grep "Estimated genome size" $err_file | awk '{print $4}' )
    coverage=$(grep "Estimated coverage" $err_file | awk '{print $3}' )
  '''
}

params.mash_reference = '/db/RefSeqSketchesDefaults.msh'
params.mash_options = '-v 0 -d 0.5'
process mash_dist {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/mash:latest'

  when:
  params.mash

  input:
  tuple val(sample), file(msh) from mash_sketch_files.concat(fastas_mash)

  output:
  tuple sample, file("mash/${sample}_mashdist.txt")
  tuple sample, env(genus) into mash_genus_results, mash_genus_prokka, mash_genus_gc, mash_genus_amrfinder
  tuple sample, env(species) into mash_species_results, mash_species_prokka, mash_species_gc, mash_species_amrfinder
  tuple sample, env(full_mash) into mash_full_results
  tuple sample, env(pvalue) into mash_pvalue_results
  tuple sample, env(distance) into mash_distance_results
  tuple sample, env(salmonella_flag) into salmonella_flag
  tuple sample, env(ecoli_flag) into ecoli_flag, shigella_flag
  tuple sample, env(klebsiella_flag) into klebsiella_flag
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p mash logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mash dist -p !{task.cpus} !{params.mash_options} !{params.mash_reference} !{msh} | sort -gk3 > mash/!{sample}_mashdist.txt 2>> $err_file

    if [ ! -s "mash/!{sample}_mashdist.txt" ]
    then
      echo "!{sample} had no mash results with '!{params.mash_options}'. Trying again without those parameters."
      mash dist -p !{task.cpus} !{params.mash_reference} !{msh} | sort -gk3 > mash/!{sample}_mashdist.txt 2>> $err_file
    fi

    mash_result=($(head -n 1 mash/!{sample}_mashdist.txt | head -n 1 | cut -f 1 | cut -f 8 -d "-" | cut -f 1,2 -d "_" | cut -f 1 -d "." | tr "_" " " ) 'missing' 'missing')
    genus=${mash_result[0]}
    species=${mash_result[1]}
    full_mash=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 1 )
    pvalue=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 4 )
    distance=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 3 )
    if [ -z "$full_mash" ] ; then full_mash='missing' ; fi
    if [ -z "$pvalue" ] ; then pvalue='NA' ; fi
    if [ -z "$distance" ] ; then distance='NA' ; fi

    find_salmonella=''
    salmonella_flag=''
    find_salmonella=$(head mash/!{sample}_mashdist.txt | grep "Salmonella" | head -n 1)
    if [ -n "$find_salmonella" ]
    then
      salmonella_flag="found"
    else
      salmonella_flag="not"
    fi

    find_ecoli=''
    ecoli_flag=''
    find_ecoli=$(head mash/!{sample}_mashdist.txt | grep -e "Escherichia" -e "Shigella" | head -n 1)
    if [ -n "$find_ecoli" ]
    then
      ecoli_flag="found"
    else
      ecoli_flag="not"
    fi

    find_klebsiella=''
    klebsiella_flag=''
    find_klebsiella=$(head mash/!{sample}_mashdist.txt | grep -e "Klebsiella" -e "Enterobacter" -e "Serratia" | head -n 1)
    if [ -n "$find_klebsiella" ]
    then
      klebsiella_flag="found"
    else
      klebsiella_flag="not"
    fi
  '''
}

params.prokka = false
params.prokka_options = ''
params.center = 'STAPHB'
params.mincontiglen = 500
process prokka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.maxcpus
  container 'staphb/prokka:latest'

  when:
  params.prokka

  input:
  tuple val(sample), file(contigs), val(genus), val(species) from contigs_prokka.concat(fastas_prokka).join(mash_genus_prokka, remainder: true, by:0).join(mash_species_prokka, remainder: true, by:0)

  output:
  file("prokka/${sample}/${sample}.{err,faa,ffn,fna,fsa,gbk,gff,log,sqn,tbl,tsv}")
  file("prokka/${sample}/${sample}.txt") into prokka_files
  file("gff/${sample}.gff")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p prokka gff logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    prokka -v >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    prokka !{params.prokka_options} \
      --cpu !{task.cpus} \
      --compliant \
      --centre !{params.center} \
      --mincontiglen !{params.mincontiglen} \
      --outdir prokka/!{sample} \
      --locustag locus_tag \
      --prefix !{sample} \
      --genus !{genus} \
      --species !{species} \
      --force !{contigs} 2>> $err_file | tee -a $log_file

    cp prokka/!{sample}/!{sample}.gff gff/!{sample}.gff
  '''
}

params.quast = true
params.quast_options = ''
process quast {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/quast:latest'

  when:
  params.quast

  input:
  tuple val(sample), file(contigs) from contigs_quast.concat(fastas_quast)

  output:
  file("quast/${sample}/*")
  file("quast/${sample}_quast_report.tsv") into quast_files
  file("quast/${sample}/{basic_stats,icarus_viewers}/*{pdf,html}")
  file("quast/${sample}/transposed_report.tsv") into quast_files_combine
  tuple sample, env(gc) into quast_gc_results
  tuple sample, env(num_contigs) into quast_contigs_results
  tuple sample, env(n50) into quast_N50_contigs_results
  tuple sample, env(length) into quast_length_results, quast_length_gc
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    quast.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    quast.py !{params.quast_options} \
      !{contigs} \
      --output-dir quast/!{sample} \
      --threads !{task.cpus} \
      2>> $err_file | tee -a $log_file

    gc=$(grep "GC (" quast/!{sample}/report.txt | awk '{print $3}' )
    num_contigs=$(grep "contigs" quast/!{sample}/report.txt | grep -v "(" | awk '{print $3}' )
    n50=$(grep "N50" quast/!{sample}/report.txt | awk '{print $2}' )
    length=$(grep "Total length" quast/!{sample}/report.txt | grep -v "(" | awk '{print $3}' )
    if [ -z "$gc" ] ; then gc='NA' ; fi
    if [ -z "$num_contigs" ] ; then num_contigs='NA' ; fi
    if [ -z "$n50" ] ; then n50='NA' ; fi
    if [ -z "$length" ] ; then length='NA' ; fi

    cp quast/!{sample}/report.tsv quast/!{sample}_quast_report.tsv
  '''
}

quast_files_combine
  .collectFile(name: "report.tsv",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/quast")

process shuffle {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/lyveset:latest'

  when:
  params.cg_pipeline

  input:
  tuple val(sample), file(reads) from clean_reads_cg

  output:
  tuple sample, file("shuffled/${sample}_shuffled.fastq.gz") into shuffled_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p shuffled logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    run_assembly_shuffleReads.pl -gz !{reads} 2>> $err_file > shuffled/!{sample}_shuffled.fastq.gz
  '''
}

shuffled_files
  .join(quast_length_gc, remainder: true, by:0)
  .join(mash_genome_size_gc, remainder: true, by:0)
  .join(mash_genus_gc, remainder: true, by:0)
  .join(mash_species_gc, remainder: true, by:0)
  .combine(genome_sizes)
  .set { for_gc }

params.cg_pipeline_options = '--qual_offset 33 --minLength 1'
process cg_pipeline {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/lyveset:latest'

  when:
  params.cg_pipeline

  input:
  tuple val(sample), file(fastq), val(quast), val(mash), val(genus), val(species), file(genome_file) from for_gc

  output:
  file("cg_pipeline/${sample}_cg_pipeline_report.txt") into cg_pipeline_files
  tuple sample, env(read_length) into cg_avrl_results
  tuple sample, env(quality) into cg_quality_results
  tuple sample, env(coverage) into cg_cov_results
  tuple sample, env(reference_genome_length) into ref_genome_length
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    genome_length=''
    if [ "!{genus}" != "null" ] && [ "!{species}" != "null" ] ; then genome_length=$(grep !{genus} !{genome_file} | grep !{species} | grep -v "#" | head -n 1 | cut -f 2 -d ":" | cut -f 1 -d "," | awk '{ print $0 "e+06" }') ; fi
    if [ -z "$genome_length" ] && [ "!{mash}" != "null" ] ; then genome_length=$(echo !{mash} | xargs printf "%.0f" ) ; fi
    if [ -z "$genome_length" ] && [ "!{quast}" != "null" ] ; then genome_length=!{quast} ; fi

    if [ -n "$genome_length" ]
    then
      run_assembly_readMetrics.pl !{fastq} \
        !{params.cg_pipeline_options} \
        --fast \
        --numcpus !{task.cpus} \
        -e $genome_length \
        2>> $err_file > cg_pipeline/!{sample}_cg_pipeline_report.txt

        read_length=$(cut -f 2 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
        quality=$(cut -f 6 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
        coverage=$(cut -f 9 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
    else
      genome_length='0'
      read_length='NA'
      quality='NA'
      coverage='NA'
      echo "Could not determine genome length of isolate, so could not run GC pipeline" | tee $log_file
    fi

    if [ -z "$read_length" ] ; then read_length='NA' ; fi
    if [ -z "$quality" ] ; then quality='NA' ; fi
    if [ -z "$coverage" ] ; then coverage='NA' ; fi
    reference_genome_length=$genome_length
  '''
}

cg_pipeline_files
  .collectFile(name: "cg_pipeline_report.txt",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/cg_pipeline")

params.seqsero2 = true
params.seqsero_options = '-m a -b mem'
process seqsero2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/seqsero2:latest'

  when:
  params.seqsero2 && flag =~ 'found'

  input:
  tuple val(sample), file(fastq_or_fasta), val(flag) from clean_reads_seqsero2.concat(fastas_seqsero2).join(salmonella_flag, by:0)

  output:
  tuple sample, env(antigenic_profile) into seqsero2_profile_results
  tuple sample, env(serotype) into seqsero2_serotype_results
  tuple sample, env(contamination) into seqsero2_contamination_results
  file("seqsero2/${sample}/*")
  file("seqsero2/${sample}/SeqSero_result.tsv") into seqsero2_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    SeqSero2_package.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fastq_check=$(echo "!{fastq_or_fasta}" | grep "fastq" | head -n 1 )
    if [ -n "$fastq_check" ] ; then type=2 ; else type=4 ; fi

    SeqSero2_package.py !{params.seqsero_options} \
      -t $type \
      -i !{fastq_or_fasta} \
      -p !{task.cpus} \
      -d seqsero2/!{sample} \
      -n !{sample} \
      2>> $err_file >> $log_file

    serotype=$(cut -f 9 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    contamination=$(cut -f 10 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    antigenic_profile=$(cut -f 8 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    enteritidis_check=$(grep "Enteritidis" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1 )
    sdf_check=$(grep "Detected Sdf" seqsero2/!{sample}/SeqSero_result.tsv | head -n 1 )

    if [ -n "$sdf_check" ] && [ -n "$enteritidis_check" ]
    then
      serotype="$serotype (Sdf+)"
    elif [ -z "$sdf_check" ] && [ -n "$enteritidis_check" ]
    then
      serotype="$serotype (Sdf-)"
    fi

    if [ -z "$serotype" ] ; then serotype='NA' ; fi
    if [ -z "$contamination" ] ; then contamination='NA' ; fi
    if [ -z "$antigenic_profile" ] ; then antigenic_profile='NA' ; fi
  '''
}

seqsero2_files
  .collectFile(name: "SeqSero_result.tsv",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/seqsero2")

params.shigatyper = true
params.shigatyper_options = ''
process shigatyper {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'andrewlangvt/shigatyper:1'

  when:
  params.shigatyper && flag =~ 'found'

  input:
  tuple val(sample), file(fastq), val(flag) from clean_reads_shigatyper.join(shigella_flag, by:0)

  output:
  file("${task.process}/${sample}_shigatyper.tsv")
  tuple sample, env(predictions) into shigatyper_predictions
  tuple sample, env(lacy_cada) into shigatyper_cadA
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    shigatyper --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    shigatyper !{params.shigatyper_options} \
      --name !{sample} \
      !{fastq} \
      2>> $err_file \
      > !{task.process}/!{sample}_shigatyper.tsv

    predictions=$(grep -v "prediction" !{task.process}/!{sample}_shigatyper.tsv | cut -f 2 | tr '\\n' ',' | sed 's/,$//g' )
    lacy_cada="$(grep -ie "lac" -ie "cad" $err_file | head -n 1)"
    if [ -z "$predictions" ] ; then predictions='none' ; fi
    if [ -z "$lacy_cada" ] ; then lacy_cada='none' ; fi
  '''
}

params.kleborate = true
params.kleborate_options = '-all'
process kleborate {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/kleborate:latest'

  when:
  params.kleborate && flag =~ 'found'

  input:
  tuple val(sample), file(contig), val(flag) from contigs_kleborate.concat(fastas_kleborate).join(klebsiella_flag, by:0)

  output:
  tuple sample, env(kleborate_score) into kleborate_score
  tuple sample, env(kleborate_mlst) into kleborate_mlst
  file("${task.process}/${sample}_results.txt") into kleborate_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    kleborate --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    kleborate !{params.kleborate_options} \
      -o !{task.process}/!{sample}_results.txt \
      -a !{contig} \
      2>> $err_file >> $log_file

    virulence_column=$(head -n 1 !{task.process}/!{sample}_results.txt | tr '\\t' '\\n' | grep -n virulence_score | cut -f 1 -d ":" )
    mlst_column=$(head -n 1 !{task.process}/!{sample}_results.txt | tr '\\t' '\\n' | grep -n ST | cut -f 1 -d ":" )
    if [ -n "$virulence_column" ] ; then kleborate_score=$(cut -f $virulence_column !{task.process}/!{sample}_results.txt | tail -n 1 ) ; fi
    if [ -n "$mlst_column" ] ; then kleborate_mlst=$(cut -f $mlst_column !{task.process}/!{sample}_results.txt | tail -n 1 ) ; fi
    if [ -z "$kleborate_score" ] ; then kleborate_score='NA' ; fi
    if [ -z "$kleborate_mlst" ] ; then kleborate_mlst='NA' ; fi
  '''
}

kleborate_files
  .collectFile(name: "kleborate_results.txt",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/kleborate")

params.serotypefinder = true
params.serotypefinder_options = ''
process serotypefinder {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/serotypefinder:latest'

  when:
  params.serotypefinder && flag =~ 'found'

  input:
  tuple val(sample), file(fasta), val(flag) from contigs_serotypefinder.concat(fastas_serotypefinder).join(ecoli_flag, by:0)

  output:
  file("${task.process}/${sample}/*")
  tuple sample, env(o_type) into serotypefinder_results_o
  tuple sample, env(h_type) into serotypefinder_results_h
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    serotypefinder.py !{params.serotypefinder_options} \
      -i !{fasta} \
      -o !{task.process}/!{sample} \
      -x \
      2>> $err_file >> $log_file

    h_type=$(cut -f 3 !{task.process}/!{sample}/results_tab.tsv | grep ^H | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    o_type=$(cut -f 3 !{task.process}/!{sample}/results_tab.tsv | grep ^O | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    if [ -z "$h_type" ] ; then h_type="none" ; fi
    if [ -z "$o_type" ] ; then o_type="none" ; fi
  '''
}

params.amrfinderplus = true
params.amrfinderplus_options = ''
process amrfinderplus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/ncbi-amrfinderplus:latest'

  when:
  params.amrfinderplus

  input:
  tuple val(sample), file(contigs), val(genus), val(species) from contigs_amrfinder.concat(fastas_amrfinder).join(mash_genus_amrfinder, by: 0).join(mash_species_amrfinder, by: 0)

  output:
  file("ncbi-AMRFinderplus/${sample}_amrfinder_plus.txt") into amrfinder_files
  tuple sample, env(amr_genes) into amr_genes
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p ncbi-AMRFinderplus logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    organism_options=(Acinetobacter_baumannii
    Campylobacter
    Enterococcus_faecalis
    Enterococcus_faecium
    Escherichia:Shigella
    Klebsiella
    Salmonella
    Staphylococcus_aureus
    Staphylococcus_pseudintermedius
    Streptococcus_agalactiae
    Streptococcus_pneumoniae
    Streptococcus_pyogenes
    Vibrio_cholerae)

    organism=$(history -p ${organism_options[@]} | grep -i !{genus} | grep -i !{species} | head -n 1 )
    if [ -z "$organism" ] ; then organism=$(history -p ${organism_options[@]} | grep -i !{genus} | head -n 1 | cut -f 1 -d ":" ) ; fi
    if [ -n "$organism" ]
    then
      organism_check="--organism $organism"
      echo "Mash result of !{genus} !{species} matched with $organism" >> $log_file
    else
      organism_check=''
      echo "Mash result of !{genus} !{species} did not match any of the organisms" >> $log_file
    fi

    amrfinder !{params.amrfinderplus_options} \
      --nucleotide !{contigs} \
      --threads !{task.cpus} \
      --name !{sample} \
      --output ncbi-AMRFinderplus/!{sample}_amrfinder_plus.txt \
      $organism_check \
      --plus \
      2>> $err_file >> $log_file

    amr_genes=$(cut -f 7 ncbi-AMRFinderplus/!{sample}_amrfinder_plus.txt | tail +2 | tr '\\n' ',' | sed 's/,$//g' )
    if [ -z "$amr_genes" ] ; then amr_genes="none" ; fi
  '''
}

amrfinder_files
  .collectFile(name: "ncbi-AMRFinderplus.txt",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/ncbi-AMRFinderplus")

params.kraken2_options = ''
process kraken2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.maxcpus
  container 'staphb/kraken2:latest'

  when:
  params.kraken2

  input:
  tuple val(sample), file(clean), path(kraken2_db) from clean_reads_kraken2.combine(local_kraken2)

  output:
  file("${task.process}/${sample}_kraken2_report.txt") into kraken2_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(top_hit) into kraken2_top_hit
  tuple sample, env(top_perc) into kraken2_top_perc
  tuple sample, env(top_reads) into kraken2_top_reads

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    kraken2 --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    kraken2 !{params.kraken2_options} \
      --paired \
      --classified-out cseqs#.fq \
      --threads !{task.cpus} \
      --db !{kraken2_db} \
      !{clean} \
      --report !{task.process}/!{sample}_kraken2_report.txt \
      2>> $err_file >> $log_file

    top_hit=$(cat !{task.process}/!{sample}_kraken2_report.txt   | grep -w S | sort | tail -n 1 | awk '{print $6 " " $7}')
    top_perc=$(cat !{task.process}/!{sample}_kraken2_report.txt  | grep -w S | sort | tail -n 1 | awk '{print $1}')
    top_reads=$(cat !{task.process}/!{sample}_kraken2_report.txt | grep -w S | sort | tail -n 1 | awk '{print $3}')
    if [ -z "$top_hit" ] ; then top_hit="NA" ; fi
    if [ -z "$top_perc" ] ; then top_perc="0" ; fi
    if [ -z "$top_reads" ] ; then top_reads="0" ; fi
  '''
}

params.local_db_type = 'nt'
process blastn {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'ncbi/blast:latest'

  when:
  params.blobtools

  input:
  tuple val(sample), file(contig), path(blastdb) from contigs_blastn.combine(local_blastdb)

  output:
  tuple sample, file("blastn/${sample}.tsv") into blastn_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    blastn -version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blastn -query !{contig} \
      -out blastn/!{sample}.tsv \
      -num_threads !{task.cpus} \
      -db !{blastdb}/!{params.local_db_type} \
      -outfmt '6 qseqid staxids bitscore std' \
      -max_target_seqs 10 \
      -max_hsps 1 \
      -evalue 1e-25 \
      2>> $err_file >> $log_file
  '''
}

params.bwa_options = ''
process bwa {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/bwa:latest'

  when:
  params.blobtools

  input:
  tuple val(sample), file(reads), file(contig) from clean_reads_bwa.join(contigs_bwa, by:0)

  output:
  tuple sample, file("${task.process}/${sample}.sam") into sam_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    bwa index !{contig} 2>> $err_file >> $log_file

    bwa mem !{params.bwa_options} \
      -t !{task.cpus} \
      !{contig} \
      !{reads} \
      2>> $err_file \
      > !{task.process}/!{sample}.sam
  '''
}

params.samtools_sort_options=''
process sort {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/samtools:latest'

  when:
  params.blobtools

  input:
  tuple val(sample), file(sam) from sam_files

  output:
  tuple sample, file("aligned/${sample}.sorted.bam*") into bam_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p aligned logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    samtools --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    samtools sort !{params.samtools_sort_options} \
      --threads !{task.cpus} \
      !{sam} \
      -o aligned/!{sample}.sorted.bam \
      --write-index \
      2>> $err_file >> $log_file
  '''
}

params.blobtools_create_options=''
process create {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'chrishah/blobtools:v1.1.1'
  when:
  params.blobtools

  input:
  tuple val(sample), file(contig), file(blastn), file(bam) from contigs_create.join(blastn_files, by:0).join(bam_files, by:0)

  output:
  tuple sample, file("blobtools/${sample}.blobDB.json") into create_files_view, create_files_plot
  file("blobtools/${sample}.${sample}.sorted.bam.cov")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools create !{params.blobtools_create_options} \
      -o blobtools/!{sample} \
      -i !{contig} \
      -b !{sample}.sorted.bam \
      -t !{blastn} \
      2>> $err_file >> $log_file
  '''
}

params.blobtools_view_options=''
process view {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'chrishah/blobtools:v1.1.1'
  when:
  params.blobtools

  input:
  tuple val(sample), file(json) from create_files_view

  output:
  tuple sample, file("blobtools/${sample}.blobDB.table.txt") into view_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools view !{params.blobtools_view_options} \
      -i !{json} \
      -o blobtools/ \
      2>> $err_file >> $log_file
  '''
}

params.blobtools_plot_options = ''
params.blobtools_resolution = 'species'
params.blobtools_format = 'png'
process blobtools {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'chrishah/blobtools:v1.1.1'

  when:
  params.blobtools

  input:
  tuple val(sample), file(json) from create_files_plot

  output:
  file("blobtools/${sample}.*.blobplot.bam0.png")
  file("blobtools/${sample}.*.blobplot.read_cov.bam0.png")
  file("blobtools/${sample}.*.blobplot.stats.txt")
  tuple sample, env(blobtools_species) into blobtools_species_results
  tuple sample, env(blobtools_perc) into blobtools_perc_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools plot !{params.blobtools_plot_options} \
      -i !{json} \
      -o blobtools/ \
      -r !{params.blobtools_resolution} \
      --format !{params.blobtools_format} \
      2>> $err_file 2>> $log_file

    perc='0.0'
    blobtools_species='missing'
    while read line
    do
      new_perc=$(echo $line | cut -f 13 -d " " | sed 's/%//g')
      min=$(echo $perc $new_perc | awk '{if ($1 > $2) print $1; else print $2}')
      if [ "$min" != "$perc" ]
      then
        perc=$new_perc
        blobtools_species=$(echo $line | cut -f 1 -d " " )
        blobtools_perc=$(echo $line | cut -f 13 -d " " )
      fi
    done < <(grep -vw all blobtools/!{sample}*.stats.txt | grep -v "# name" | tr ' ' '_' | grep '%')
  '''
}

params.mlst = true
process mlst {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/mlst:latest'

  when:
  params.mlst

  input:
  tuple val(sample), file(contig) from contigs_mlst.concat(fastas_mlst)

  output:
  file("${task.process}/${sample}_mlst.txt") into mlst_files
  tuple sample, env(mlst) into mlst_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    mlst --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mlst !{contig} 2>> $err_file > mlst/!{sample}_mlst.txt

    mlst=$(awk '{ print $2 ":" $3 }' mlst/!{sample}_mlst.txt)
  '''
}

mlst_files
  .collectFile(name: "mlst_result.tsv",
    keepHeader: false,
    sort: true,
    storeDir: "${params.outdir}/mlst")

reads
  .concat(fastas)
  .join(seqyclean_perc_kept_results, remainder: true, by: 0)
  .join(seqyclean_pairskept_results, remainder: true, by: 0)
  .join(fastqc_1_results, remainder: true, by: 0)
  .join(fastqc_2_results, remainder: true, by: 0)
  .join(mash_genome_size_results, remainder: true, by: 0)
  .join(mash_coverage_results, remainder: true, by: 0)
  .join(mash_genus_results, remainder: true, by: 0)
  .join(mash_species_results, remainder: true, by: 0)
  .join(mash_full_results, remainder: true, by: 0)
  .join(mash_pvalue_results, remainder: true, by: 0)
  .join(mash_distance_results, remainder: true, by: 0)
  .join(quast_gc_results, remainder: true, by: 0)
  .join(quast_contigs_results, remainder: true, by: 0)
  .join(quast_N50_contigs_results, remainder: true, by: 0)
  .join(quast_length_results, remainder: true, by: 0)
  .join(cg_avrl_results, remainder: true, by: 0)
  .join(cg_quality_results, remainder: true, by: 0)
  .join(cg_cov_results, remainder: true, by: 0)
  .join(ref_genome_length, remainder: true, by: 0)
  .join(seqsero2_profile_results, remainder: true, by: 0)
  .join(seqsero2_serotype_results, remainder: true, by: 0)
  .join(seqsero2_contamination_results, remainder: true, by: 0)
  .join(serotypefinder_results_o, remainder: true, by: 0)
  .join(serotypefinder_results_h, remainder: true, by: 0)
  .join(shigatyper_predictions, remainder: true, by: 0)
  .join(shigatyper_cadA, remainder: true, by: 0)
  .join(kleborate_score, remainder: true, by: 0)
  .join(kleborate_mlst, remainder: true, by: 0)
  .join(blobtools_species_results, remainder: true, by: 0)
  .join(blobtools_perc_results, remainder: true, by: 0)
  .join(kraken2_top_hit, remainder: true, by: 0)
  .join(kraken2_top_perc, remainder: true, by: 0)
  .join(kraken2_top_reads, remainder: true, by: 0)
  .join(amr_genes, remainder: true, by: 0)
  .join(mlst_results, remainder: true, by: 0)
  .set { results }

process summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${sample}"
  cpus 1
  container 'staphb/parallel-perl:latest'

  input:
  set val(sample), file(file),
    val(seqyclean_perc_kept_results),
    val(seqyclean_pairskept_results),
    val(fastqc_1_results),
    val(fastqc_2_results),
    val(mash_genome_size_results),
    val(mash_coverage_results),
    val(mash_genus_results),
    val(mash_species_results),
    val(mash_full_results),
    val(mash_pvalue_results),
    val(mash_distance_results),
    val(quast_gc_results),
    val(quast_contigs_results),
    val(quast_N50_contigs_results),
    val(quast_length_results),
    val(cg_avrl_results),
    val(cg_quality_results),
    val(cg_cov_results),
    val(ref_genome_length),
    val(seqsero2_profile_results),
    val(seqsero2_serotype_results),
    val(seqsero2_contamination_results),
    val(serotypefinder_results_o),
    val(serotypefinder_results_h),
    val(shigatyper_predictions),
    val(shigatyper_cadA),
    val(kleborate_score),
    val(kleborate_mlst),
    val(blobtools_species_results),
    val(blobtools_perc_results),
    val(kraken2_top_hit),
    val(kraken2_top_perc),
    val(kraken2_top_reads),
    val(amr_genes),
    val(mlst_results) from results

  output:
  file("summary/${sample}.summary.txt") into summary_files_txt
  file("summary/${sample}.summary.tsv") into summary_files_tsv
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    sample_id_split=($(echo !{sample} | sed 's/-/ /g' | sed 's/_/ /g' ))
    if [ "${#sample_id_split[@]}" == "5" ]
    then
      sample_id="${sample_id_split[0]}-${sample_id_split[1]}"
    else
      sample_id=${sample_id_split[0]}
    fi

    header="sample_id;sample"
    result="$sample_id;!{sample}"

    header=$header";seqyclean_pairs_kept;seqyclean_percent_kept"
    result=$result";!{seqyclean_pairskept_results};!{seqyclean_perc_kept_results}"

    header=$header";fastqc_1_reads;fastqc_2_reads"
    result=$result";!{fastqc_1_results};!{fastqc_2_results}"

    header=$header";mash_genome_size;mash_coverage;mash_genus;mash_species;mash_full;mash_pvalue;mash_distance"
    result=$result";!{mash_genome_size_results};!{mash_coverage_results};!{mash_genus_results};!{mash_species_results};!{mash_full_results};!{mash_pvalue_results};!{mash_distance_results}"

    header=$header";quast_gc_%;quast_contigs;quast_N50;quast_length"
    result=$result";!{quast_gc_results};!{quast_contigs_results};!{quast_N50_contigs_results};!{quast_length_results}"

    header=$header";cg_average_read_length;cg_average_quality;cg_coverage;ref_genome_length"
    result=$result";!{cg_avrl_results};!{cg_quality_results};!{cg_cov_results};!{ref_genome_length}"

    header=$header";seqsero2_profile;seqsero2_serotype;seqsero2_contamination"
    result=$result";!{seqsero2_profile_results};!{seqsero2_serotype_results};!{seqsero2_contamination_results}"

    header=$header";serotypefinder_o_group;serotypefinder_h_group"
    result=$result";!{serotypefinder_results_o};!{serotypefinder_results_h}"

    header=$header";shigatyper_predictions;shigatyper_cadA"
    result=$result";!{shigatyper_predictions};!{shigatyper_cadA}"

    header=$header";kleborate_score;kleborate_mlst"
    result=$result";!{kleborate_score};!{kleborate_mlst}"

    header=$header";amr_genes"
    result=$result";!{amr_genes}"

    header=$header";blobtools_top_species;blobtools_percentage"
    result=$result";!{blobtools_species_results};!{blobtools_perc_results}"

    header=$header";kraken2_top_species;kraken2_num_reads;kraken2_percentage"
    result=$result";!{kraken2_top_hit};!{kraken2_top_reads};!{kraken2_top_perc}"

    header=$header";mlst"
    result=$result";!{mlst_results}"

    echo $header > summary/!{sample}.summary.txt
    echo $result >> summary/!{sample}.summary.txt

    cat summary/!{sample}.summary.txt | tr ';' '\t' > summary/!{sample}.summary.tsv
  '''
}

summary_files_txt
  .collectFile(name: "grandeur_summary.txt",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/summary")

summary_files_tsv
  .collectFile(name: "grandeur_results.tsv",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}")
  .subscribe { println("Summary can be found at $it") }

params.multiqc_options = ''
params.multiqc = true
process multiqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "multiqc"
  cpus 1
  container 'ewels/multiqc:latest'

  when:
  params.multiqc

  input:
  file(fastqc) from fastqc_files.collect().ifEmpty([])
  file(quast) from quast_files.collect().ifEmpty([])
  file(seqyclean) from seqyclean_files.collect().ifEmpty([])
  file(kraken2) from kraken2_files.concat(kraken2_files_empty).collect().ifEmpty([])
  file(prokka) from prokka_files.collect().ifEmpty([])

  output:
  file("${task.process}/multiqc_report.html")
  file("${task.process}/multiqc_data/*")
  file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} quast fastqc kraken2 prokka seqyclean logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    multiqc --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    for quast_file in !{quast}
    do
      sample=$(echo $quast_file | sed 's/_quast_report.tsv//g' | head -n 1 )
      mkdir -p quast/$sample
      mv $quast_file quast/$sample/report.tsv
    done

    multiqc !{params.multiqc_options} \
      --outdir !{task.process} \
      --cl_config "prokka_fn_snames: True"  \
      . \
      2>> $err_file >> $log_file
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("MultiQC report can be found at $params.outdir/multiqc/multiqc_report.html")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}

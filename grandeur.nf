#!/usr/bin/env nextflow

println("Currently using the Grandeur workflow for use with microbial Illumina MiSeq sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210430")
println("")

//TBA : bgzip the seqyclean files?

params.reads = workflow.launchDir + '/reads'
params.outdir = workflow.launchDir + '/grandeur'

params.maxmem = Math.round(Runtime.runtime.totalMemory() / 10241024).GB
println("The maximum amount of memory used in this workflow is ${params.maxmem}")
params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq*",
                  "${params.reads}/*_{1,2}.fastq*"], size: 2 )
  .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .ifEmpty{
    println("FATAL : No fastq or fastq.gz files were found at ${params.reads}")
    println("Set 'params.reads' to directory with paired-end reads")
    exit 1
  }
  .into { reads_seqyclean ; reads_fastqc }

params.cg_pipeline = 'true'
params.genome_sizes = workflow.projectDir + "/configs/genome_sizes.json"
genome_sizes = params.cg_pipeline
              ? Channel.fromPath(params.genome_sizes, type:'file', checkIfExists: true).view { "Genome Sizes File : $it" }
              : Channel.empty()

params.blobtools = false
params.local_blastdb='/blast/blastdb'
local_blastdb = params.blobtools
              ? Channel.fromPath(params.local_blastdb, type:'dir').view { "Local Blast Database : $it" }
              : Channel.empty()

params.seqyclean_minlen = 25
params.seqyclean_contaminant_file = "/Adapters_plus_PhiX_174.fasta"
params.seqyclean_options = ''
process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  container 'staphb/seqyclean:latest'

  input:
  set val(sample), file(reads) from reads_seqyclean

  output:
  tuple sample, file("seqyclean/${sample}_clean_PE{1,2}.fastq.gz") into clean_reads_shovill, clean_reads_mash, clean_reads_cg, clean_reads_seqsero2, clean_reads_bwa, clean_reads_shigatyper, clean_reads_serotypefinder
  file("seqyclean/${sample}_clean_SE.fastq.gz")
  file("seqyclean/${sample}_clean_SummaryStatistics.{txt,tsv}")
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
    echo "seqyclean version: $(seqyclean -h | grep Version | head -n 1)" >> $log_file

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

params.shovill_options = ''
process shovill {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  memory params.maxmem
  cpus params.maxcpus
  errorStrategy { task.attempt<=3 ? 'retry' : 'ignore' }
  container 'staphb/shovill:latest'

  input:
  set val(sample), file(reads) from clean_reads_shovill

  output:
  file("shovill/${sample}/contigs.{fa,gfa}")
  file("shovill/${sample}/shovill.{corrections,log}")
  file("shovill/${sample}/spades.fasta")
  tuple sample, file("contigs/${sample}_contigs.fa") into contigs_prokka, contigs_quast, contigs_gc, contigs_abricate, contigs_abricate_ecoli, contigs_index, contigs_blastn, contigs_mlst, contigs_bwa, contigs_create, contigs_kleborate, contigs_amrfinder
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} contigs logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    shovill --version >> $log_file

    echo "the cpus are !{task.cpus}"
    echo "the memory is !{task.memory}"

    shovill !{params.shovill_options}\
      --cpu !{task.cpus} \
      --ram !{task.memory} \
      --outdir !{task.process}/!{sample} \
      --R1 !{reads[0]} \
      --R2 !{reads[1]} \
      2>> $err_file >> tee $log_file
    cp shovill/!{sample}/contigs.fa contigs/!{sample}_contigs.fa
  '''
}

params.fastqc = true
params.fastqc_options = ''
process fastqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  container 'staphb/fastqc:latest'

  when:
  params.fastqc

  input:
  set val(sample), file(raw) from reads_fastqc

  output:
  file("fastqc/*.{html,zip}")
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
    fastqc --version >> $log_file

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
  echo false
  cpus 1
  container 'staphb/mash:latest'

  when:
  params.mash

  input:
  set val(sample), file(reads) from clean_reads_mash

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
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file

    cat !{reads} | mash sketch -m 2 -o mash/!{sample} - 2>> $err_file | tee $log_file

    genome_size=$(grep "Estimated genome size" $err_file | awk '{print $4}')
    coverage=$(grep "Estimated coverage" $err_file | awk '{print $3}')
  '''
}

params.mash_reference = "/db/RefSeqSketchesDefaults.msh"
params.mash_pvalue = 0
params.mash_distance = 0.5
params.mash_options = ''
process mash_dist {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/mash:latest'

  when:
  params.mash

  input:
  set val(sample), file(msh) from mash_sketch_files

  output:
  tuple sample, file("mash/${sample}_mashdist.txt") into mash_dist_files
  tuple sample, env(genus) into mash_genus_results, mash_genus_prokka, mash_genus_gc
  tuple sample, env(species) into mash_species_results, mash_species_prokka, mash_species_gc
  tuple sample, env(full_mash) into mash_full_results
  tuple sample, env(pvalue) into mash_pvalue_results
  tuple sample, env(distance) into mash_distance_results
  tuple sample, env(salmonella_flag) into salmonella_flag
  tuple sample, env(ecoli_flag) into ecoli_flag, ecoli_flag2, shigella_flag
  tuple sample, env(klebsiella_flag) into klebsiella_flag
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p mash logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file

    mash dist -p !{task.cpus} !{params.mash_options} -v !{params.mash_pvalue} -d !{params.mash_distance} !{params.mash_reference} !{msh} | sort -gk3 > mash/!{sample}_mashdist.txt 2>> $err_file

    if [ ! -s "mash/!{sample}_mashdist.txt" ]
    then
      echo "!{sample} had no mash results with '-v ${params.mash_pvalue}' and -d '{params.mash_distance}'. Trying again without those parameters."
      mash dist -p !{task.cpus} !{params.mash_options} !{params.mash_reference} !{msh} | sort -gk3 > mash/!{sample}_mashdist.txt 2>> $err_file
    fi

    mash_result=($(head -n 1 mash/!{sample}_mashdist.txt | head -n 1 | cut -f 1 | cut -f 8 -d "-" | cut -f 1,2 -d "_" | cut -f 1 -d "." | tr "_" " " ) 'missing' 'missing')
    genus=${mash_result[0]}
    species=${mash_result[1]}
    full_mash=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 1)
    pvalue=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 4)
    distance=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 3)

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

contigs_prokka
  .join(mash_genus_prokka, remainder: true, by:0)
  .join(mash_species_prokka, remainder: true, by:0)
  .set { for_prokka }

params.prokka = true
params.prokka_options = ''
params.center = 'STAPHB'
params.mincontiglen = 500
process prokka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  container 'staphb/prokka:latest'

  when:
  params.prokka

  input:
  set val(sample), file(contigs), val(genus), val(species) from for_prokka

  output:
  file("prokka/${sample}/${sample}.{err,faa,ffn,fna,fsa,gbk,gff,log,sqn,tbl,tsv,txt}")
  file("gff/${sample}.gff")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p prokka gff logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    prokka -v >> $log_file

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
  echo false
  cpus 1
  container 'staphb/quast:latest'

  when:
  params.quast

  input:
  set val(sample), file(contigs) from contigs_quast

  output:
  file("quast/${sample}/*")
  file("quast/${sample}/{basic_stats,icarus_viewers}/*{pdf,html}")
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
    quast.py --version >> $log_file

    quast.py !{params.quast_options}\
      !{contigs} \
      --output-dir quast/!{sample} \
      --threads !{task.cpus} 2>> $err_file | tee -a $log_file

    gc=$(grep "GC (" quast/!{sample}/report.txt | awk '{print $3}' )
    num_contigs=$(grep "contigs" quast/!{sample}/report.txt | grep -v "(" | awk '{print $3}' )
    n50=$(grep "N50" quast/!{sample}/report.txt | awk '{print $2}' )
    length=$(grep "Total length" quast/!{sample}/report.txt | grep -v "(" | awk '{print $3}' )
  '''
}

process shuffle {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  container 'staphb/lyveset:latest'

  when:
  params.cg_pipeline

  input:
  set val(sample), file(reads) from clean_reads_cg

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

params.cg_pipeline_qual_offset = 33
params.cg_pipeline_minlength = 1
params.cg_pipeline_options = ''
process cg_pipeline {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/lyveset:latest'

  when:
  params.cg_pipeline

  input:
  set val(sample), file(fastq), val(quast), val(mash), val(genus), val(species), file(genome_file) from for_gc

  output:
  file("cg_pipeline/${sample}_cg_pipeline_report.txt")
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

    genome_length=''
    if [ "!{genus}" != "null" ] && [ "!{species}" != "null" ]
    then
      genome_length=$(grep !{genus} !{genome_file} | grep !{species} | grep -v "#" | head -n 1 | cut -f 2 -d ":" | cut -f 1 -d "," | awk '{ print $0 "e+06" }')
      reference_genome_length=$genome_length
    fi

    if [ -z "$genome_length" ] && [ "!{mash}" != "null" ] ; then genome_length=$(echo !{mash} | xargs printf "%.0f" ) ; fi
    if [ -z "$genome_length" ] && [ "!{quast}" != "null" ] ; then genome_length=!{quast} ; fi

    run_assembly_readMetrics.pl !{fastq} \
      !{params.cg_pipeline_options} \
      --fast \
      --numcpus !{task.cpus} \
      -e $genome_length \
      --qual_offset !{params.cg_pipeline_qual_offset} \
      --minLength !{params.cg_pipeline_minlength} \
      2>> $err_file > cg_pipeline/!{sample}_cg_pipeline_report.txt

    read_length=$(cut -f 2 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
    quality=$(cut -f 6 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
    coverage=$(cut -f 9 cg_pipeline/!{sample}_cg_pipeline_report.txt | tail -n 1 )
  '''
}

clean_reads_seqsero2
  .join(salmonella_flag, by:0)
  .set { for_seqsero2 }

params.seqsero2 = true
params.seqsero_workflow = 'a'
params.seqsero_mapping_algorithm = 'mem'
params.seqsero_options = ''
process seqsero2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/seqsero2:latest'

  when:
  params.seqsero2 && flag =~ 'found'

  input:
  set val(sample), file(fastq), val(flag) from for_seqsero2

  output:
  tuple sample, env(antigenic_profile) into seqsero2_profile_results
  tuple sample, env(serotype) into seqsero2_serotype_results
  tuple sample, env(contamination) into seqsero2_contamination_results
  file("seqsero2/${sample}/*")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    SeqSero2_package.py --version >> $log_file

    SeqSero2_package.py !{params.seqsero_options} \
      -m !{params.seqsero_workflow} \
      -t 2 \
      -i !{fastq} \
      -p !{task.cpus} \
      -b !{params.seqsero_mapping_algorithm} \
      -d seqsero2/!{sample} \
      -n !{sample} \
      2>> $err_file >> $log_file

    antigenic_profile=$(cut -f 8 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    serotype=$(cut -f 9 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
    contamination=$(cut -f 10 seqsero2/!{sample}/SeqSero_result.tsv | tail -n 1)
  '''
}

clean_reads_shigatyper
  .join(shigella_flag, by:0)
  .set { for_shigatyper }

params.shigatyper = false
params.shigatyper_options = ''
process shigatyper {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/shigatyper:latest'

  when:
  params.shigatyper && flag =~ 'found'

  input:
  set val(sample), file(fastq), val(flag) from for_shigatyper

  output:
  file("${task.process}/${sample}/*")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    shigatyper.py --version >> $log_file

    shigatyper.py !{params.shigatyper_options} \
      -n !{sample} \
      !{fastq} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

contigs_kleborate
  .join(klebsiella_flag, by:0)
  .set { for_kleborate }

params.kleborate = false
params.kleborate_options = '-all'
process kleborate {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/kleborate:latest'

  when:
  params.kleborate && flag =~ 'found'

  input:
  set val(sample), file(contig), val(flag) from for_kleborate

  output:
  tuple sample, env(kleborate_score) into kleborate_results
  file("${task.process}/${sample}_results.txt")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    kleborate --version >> $log_file

    kleborate !{params.kleborate_options} \
      -o !{task.process}/!{sample}_results.txt \
      -a !{contig} \
      2>> $err_file >> $log_file

    kleborate_score=''
    exit 1
  '''
}

params.abricate = true
params.abricate_db = ['ncbi', 'vfdb']
params.abricate_minid = 80
params.abricate_mincov = 80
process abricate {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample} : ${db}"
  echo false
  cpus params.medcpus
  container 'staphb/abricate:latest'

  when:
  params.abricate

  input:
  set val(sample), file(contig) from contigs_abricate
  each db from params.abricate_db

  output:
  tuple sample, file("abricate/${db}/${sample}_${db}_results.txt") into abricate_files
  file("logs/${task.process}/${sample}_${db}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{db} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}_!{db}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}_!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    abricate --version >> $log_file
    abricate --list >> $log_file

    abricate \
      --db !{db} \
      --threads !{task.cpus} \
      --minid !{params.abricate_minid} \
      --mincov !{params.abricate_mincov} \
      !{contig} > abricate/!{db}/!{sample}_!{db}_results.txt \
      2>> $err_file
  '''
}

abricate_files
  .groupTuple()
  .set { abricate_results }

contigs_abricate_ecoli
  .join(ecoli_flag2, by:0)
  .set { for_abricate_ecoli }

clean_reads_serotypefinder
  .join(ecoli_flag, by:0)
  .set { for_serotypefinder }

params.serotypefinder = false
params.serotypefinder_options = ''
process serotypefinder {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/serotypefinder:latest'

  when:
  params.serotypefinder && flag =~ 'found'

  input:
  set val(sample), file(fastq), val(flag) from for_serotypefinder

  output:
  file("${task.process}/${sample}/*")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    serotypefinder --version >> $log_file

    serotypefinder !{params.serotypefinder} \
      --infile !{fastq} \
      --outputPath !{task.process} \
      --mincov !{params.serotypefinder_mincov} \
      --threshold !{params.serotypefinder_threshold} \
      --extented_output \
      2>> $err_file >> $log_file

      exit 1
  '''
}




params.abricate_ecoli = true
params.abricate_ecoli_db = ['ecoh', 'serotypefinder']
process abricate_ecoli_serotyping {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample} : ${db}"
  echo false
  cpus params.medcpus
  container 'staphb/abricate:0.8.13s'

  when:
  params.abricate_ecoli && flag =~ 'found'

  input:
  set val(sample), file(contig), val(flag) from for_abricate_ecoli
  each db from params.abricate_ecoli_db

  output:
  tuple sample, file("abricate/${db}/${sample}_${db}_results.txt") into abricate_ecoli_files
  file("logs/${task.process}/${sample}_${db}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p abricate/!{db} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}_!{db}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}_!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    abricate --version >> $log_file
    abricate --list >> $log_file

    abricate \
      --db !{db} \
      --threads !{task.cpus} \
      --minid !{params.abricate_minid} \
      --mincov !{params.abricate_mincov} \
      !{contig} > abricate/!{db}/!{sample}_!{db}_results.txt \
      2>> $err_file
  '''
}

abricate_ecoli_files
  .groupTuple()
  .set {abricate_ecoli_results}

// abricate_files.collect()
//   .concat(abricate_ecoli_files)
//   .set { for_abricate_summary }
//
// if (params.abricate && params.abricate_ecoli) {
//   abricate_summary_db = params.abricate_ecoli_db + params.abricate_db
//   } else if (params.abricate) {
//     abricate_summary_db = params.abricate_db
//   } else if (params.abricate_ecoli) {
//     abricate_summary_db = params.abricate_ecoli_db
//   } else {
//     abricate_summary_db = ''
//   }
//
// process abricate_summary {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "${db}"
//   echo false
//   cpus params.medcpus
//   container 'staphb/abricate:latest'
//
//   input:
//   file(files) from for_abricate_summary.collect()
//   each db from abricate_summary_db
//
//   output:
//   file("abricate/${db}/${db}_summary.txt")
//   file("logs/${task.process}/${db}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     mkdir -p abricate/!{db} logs/!{task.process}
//     log_file=logs/!{task.process}/!{db}.!{workflow.sessionId}.log
//     err_file=logs/!{task.process}/!{db}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     abricate --version >> $log_file
//     abricate --list >> $log_file
//
//     num=$(ls *_results.txt | grep !{db} | wc -l )
//     if [ -z "$num" ] ; then num=0 ; fi
//     if [[ "$num" > 0 ]]
//     then
//       abricate --summary *_!{db}_results.txt > abricate/!{db}/!{db}_summary.txt 2>> $err_file
//       cat *_!{db}_results.txt > abricate/!{db}/!{db}_all.txt 2>> $err_file
//     fi
//   '''
// }

params.amrfinderplus = false
params.amrfinderplus_options = ''
process amrfinderplus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample} : ${db}"
  echo true
  cpus params.medcpus
  container 'staphb/ncbi-amrfinderplus:latest'

  when:
  params.amrfinderplus

  input:
  set val(sample), file(contigs) from contigs_amrfinder

  output:
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p ncbi-AMRFinderplus/ logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}_!{db}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}_!{db}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    amrfinder -l >> $log_file

    amrfinder -h

    amrfinder -l

    amrfinder !{params.amrfinderplus_options} \
      -n !{contigs} \
      --threads !{task.cpus} \
      --name !{sample} \
      -o whatever.txt \
      --oranism tbd \
      2>> $err_file

    exit 1
  '''
}

contigs_blastn
  .combine(local_blastdb)
  .set { for_blastn }

params.local_db_type = 'nt'
process blastn {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'ncbi/blast:latest'

  when:
  params.blobtools

  input:
  set val(sample), file(contig), path(blastdb) from for_blastn

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
    blastn -version >> $log_file

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

clean_reads_bwa
  .join(contigs_bwa, by:0)
  .set { for_bwa }

params.bwa_options = ''
process bwa {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/bwa:latest'

  when:
  params.blobtools

  input:
  set val(sample), file(reads), file(contig) from for_bwa

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
    echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file

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
  echo false
  cpus params.medcpus
  container 'staphb/samtools:latest'

  when:
  params.blobtools

  input:
  set val(sample), file(sam) from sam_files

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
    samtools --version >> $log_file

    samtools sort !{params.samtools_sort_options} \
      --threads !{task.cpus} \
      !{sam} \
      -o aligned/!{sample}.sorted.bam \
      --write-index \
      2>> $err_file >> $log_file
  '''
}

contigs_create
  .join(blastn_files, by:0)
  .join(bam_files, by:0)
  .set { for_create }

params.blobtools_create_options=''
process create {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  container 'chrishah/blobtools:v1.1.1'
  when:
  params.blobtools

  input:
  set val(sample), file(contig), file(blastn), file(bam) from for_create

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
    echo "blobtools version $(blobtools -v)" >> $log_file

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
  echo false
  cpus 1
  container 'chrishah/blobtools:v1.1.1'
  when:
  params.blobtools

  input:
  set val(sample), file(json) from create_files_view

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
    echo "blobtools version $(blobtools -v)" >> $log_file

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
  echo false
  cpus 1
  container 'chrishah/blobtools:v1.1.1'

  when:
  params.blobtools

  input:
  set val(sample), file(json) from create_files_plot

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
    echo "blobtools version $(blobtools -v)" >> $log_file

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
  echo false
  cpus 1
  container 'staphb/mlst:latest'

  when:
  params.mlst

  input:
  set val(sample), file(contig) from contigs_mlst

  output:
  file("${task.process}/${sample}_mlst.txt")
  tuple sample, env(mlst) into mlst_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    mlst --version >> $log_file

    mlst !{contig} 2>> $err_file > mlst/!{sample}_mlst.txt

    mlst=$(awk '{ print $2 ":" $3 }' mlst/!{sample}_mlst.txt)
  '''
}

// // TODO : multiqc
// params.multiqc = false
// params.multiqc_options = ''
// //container 'staphb/multiqc:latest'
// // TODO : amrfinder_plus
// params.amrfinder_plus = false
// params.amrfinder_plus_options = ''
// //container 'staphb/ncbi-amrfinderplus:latest'
// // TODO : frp_plasmid
// params.rfp_plasmid = false
// params.rfp_plasmid_options = ''
// //container 'none'
// // TODO : pointfinder
// params.pointfinder = false
// params.pointfinder_options = ''
// //container 'none'
// // TODO : mubsuite
// params.mobsuite = false
// params.mobsuite_options = ''
// //container 'none'
// params.sistr = false
// params.sistr_options = ''
// //container 'staphb/sistr:latest'
// params.serotypefinder = false
// params.serotypefinder_options = ''
// //container 'staphb/serotypefinder:latest'
// params.plasmidseeker = false
// params.plasmidseeker_options = ''
// //container 'staphb/plasmidseeker:latest'
// // TODO : kleborate


seqyclean_perc_kept_results
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
  .join(blobtools_species_results, remainder: true, by: 0)
  .join(blobtools_perc_results, remainder: true, by: 0)
  .join(mlst_results, remainder: true, by: 0)
  .join(abricate_results, remainder: true, by: 0)
  .join(abricate_ecoli_results, remainder: true, by: 0)
  .set { results }

process summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${sample}"
  echo false
  cpus 1
  container 'staphb/parallel-perl:latest'

  input:
  set val(sample), val(seqyclean_perc_kept_results),
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
    val(blobtools_species_results),
    val(blobtools_perc_results),
    val(mlst_results),
    file(abricate_results),
    file(abricate_serotype_results) from results

  output:
  file("summary/${sample}.summary.txt") into summary_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    sample_id=$(echo !{sample} | cut -f 1 -d "-" )

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

    header=$header";blobtools_top_species;blobtools_percentage"
    result=$result";!{blobtools_species_results};!{blobtools_perc_results}"

    header=$header";mlst"
    result=$result";!{mlst_results}"

    abricate_result=";"
    abricate_result_header=";"
    abricate_db=($(echo "!{params.abricate_db}" | tr -d '[](){},'))
    for db in ${abricate_db[@]}
    do
      found=''
      if [ -f "!{sample}_${db}_results.txt" ]
      then
        gene_col=''
        gene_col=$(head -n 1 !{sample}_${db}_results.txt | tr '\t' '\n' | grep -n "GENE" | cut -f 1 -d ":")
        if [ -z "$gene_col" ] ; then gene_col=5 ; fi

        found="$(cat !{sample}_${db}_results.txt | cut -f $gene_col | grep -v "GENE" | tr '\\n' ' ')"
      fi
      if [ -z "$found" ] ; then found='none' ; fi
      abricate_result="$abricate_result;$found"
      abricate_result_header="$abricate_result_header;$db"
    done
    abricate_result="$(echo $abricate_result | sed 's/;;//g' | sed 's/  / /g' )"
    abricate_result_header="$(echo $abricate_result_header | sed 's/;;//g')"

    serotype_result=";"
    serotype_result_header=";"
    abricate_serotype_db=($(echo "!{params.abricate_ecoli_db}" | tr -d '[](){},'))
    for db in ${abricate_serotype_db[@]}
    do
      O=''
      H=''
      if [ -f "!{sample}_${db}_results.txt" ]
      then
        gene_col=''
        gene_col=$(head -n 1 !{sample}_${db}_results.txt | tr '\t' '\n' | grep -n "GENE" | cut -f 1 -d ":")
        if [ -z "$gene_col" ] ; then gene_col=5 ; fi

        O_H=($(cat !{sample}_${db}_results.txt | cut -f $gene_col | grep -v "GENE" | tr '/' '_'))
        for gene in ${O_H[@]}
        do
          h=($(echo "$gene" | tr "-" "\\n" | tr "_" "\\n" | grep ^"H" | tr "\\n" " " | head) ' ')
          o=($(echo "$gene" | tr "-" "\\n" | tr "_" "\\n" | grep ^"O" | tr "\\n" " " | head) ' ')
          O="${O[@]}""${o[@]}"
          H="${H[@]}""${h[@]}"
        done
      fi
      if [ -z "$O" ] ; then O='none' ; fi
      if [ -z "$H" ] ; then H='none' ; fi
      O="$(history -p ${O[@]} | sort | uniq -c | tr '\\n' ',')"
      H="$(history -p ${H[@]} | sort | uniq -c | tr '\\n' ',')"
      serotype_result="$serotype_result;$O,;$H,"
      serotype_result_header="$serotype_result_header;$db : O ; $db : H"
    done
    serotype_result="$(echo $serotype_result | sed 's/;;//g' | sed 's/,,//g' | sed 's/  / /g' )"
    serotype_result_header="$(echo $serotype_result_header | sed 's/;;//g')"

    header=$header";$abricate_result_header;$serotype_result_header"
    result=$result";$abricate_result;$serotype_result"

    echo $header > summary/!{sample}.summary.txt
    echo $result >> summary/!{sample}.summary.txt
  '''
}

process combined_summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "summary.txt"
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "logs/summary/*.{log,err}"
  publishDir "${workflow.launchDir}", mode: 'copy', overwrite: true, pattern: "grandeur_results.txt"
  tag "summary"
  echo false
  cpus 1
  container 'staphb/parallel-perl:latest'

  input:
  file(summary) from summary_files.collect()

  output:
  file("summary.txt")
  file("grandeur_results.txt")
  file("logs/${task.process}/summary.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/summary.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/summary.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    cat *.summary.txt | grep "sample_id" | head -n 1 > summary.txt
    cat *.summary.txt | grep -v "sample_id" | sort | uniq >> summary.txt 2>> $err_file

    cat summary.txt | tr ';' '\t' > grandeur_results.txt
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}

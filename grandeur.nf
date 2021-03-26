#!/usr/bin/env nextflow

println("Currently using the Grandeur workflow for use with microbial Illumina MiSeq sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210430")
println("")

//TBA : bgzip the seqyclean files?

params.reads = workflow.launchDir + '/Sequencing_reads/Raw'
params.outdir = workflow.launchDir + '/grandeur'

params.fastqc = true
params.mash = true

params.mincontiglen = 500

Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq*",
                  "${params.reads}/*_{1,2}.fastq*"], size: 2 )
  .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .ifEmpty{
    println("No fastq or fastq.gz files were found at ${params.reads}")
    println("Set 'params.reads' to directory with paired-end reads")
    exit 1
  .set { reads_seqyclean ; reads_fastqc }

process seqyclean {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1

  input:
  set val(sample), file(reads) from reads_seqyclean

  output:
  tuple sample, file("seqyclean/${sample}_clean_PE{1,2}.fastq") into clean_reads_shovill, clean_reads_mash
  file("seqyclean/${sample}_clean_SE.fastq")
  file("seqyclean/${sample}_clean_SummaryStatistics.{txt,tsv}")
  file("logs/seqyclean/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(perc_kept) into seqyclean_perc_kept_results
  tuple sample, env(kept) into seqyclean_pairskept_results

  shell:
  '''
    mkdir -p seqyclean logs/seqyclean
    log_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.log
    err_file=logs/seqyclean/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file

    kept=''
    perc_kept=''

    seqyclean -minlen !{params.seqyclean_minlen} -qual -c !{params.seqyclean_contaminant_file} -1 !{reads[0]} -2 !{reads[1]} -o seqyclean/!{sample}_clean 2>> $err_file >> $log_file
    kept=$(cut -f 58 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
    perc_kept=$(cut -f 59 seqyclean/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)

    if [ -z "$kept" ] ; then kept="0" ; fi
    if [ -z "$perc_kept" ] ; then perc_kept="0" ; fi
  '''
}

process fastqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  when:
  params.fastqc

  input:
  set val(sample), file(raw) from reads_fastqc

  output:
  file("fastqc/*.{html,zip}")
  tuple sample, env(raw_1) into fastqc_1_results
  tuple sample, env(raw_2) into fastqc_2_results
  file("logs/fastqc/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p fastqc logs/fastqc
    log_file=logs/fastqc/!{sample}.!{workflow.sessionId}.log
    err_file=logs/fastqc/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastqc --version >> $log_file

    fastqc --outdir fastqc --threads !{task.cpus} !{raw} 2>> $err_file >> $log_file

    zipped_fastq=($(ls fastqc/*fastqc.zip) "")

    raw_1=$(unzip -p ${zipped_fastq[0]} */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=$(unzip -p fastqc/*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )

    if [ -z "$raw_1" ] ; then raw_1="0" ; fi
    if [ -z "$raw_2" ] ; then raw_2="0" ; fi
  '''
}

process shovill {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  memory params.maxmem
  cpus params.maxcpus
  errorStrategy 'retry' 3 times

  input:
  set val(sample), file(reads) from clean_reads_shovill

  output:
  file("shovill/${sample}/contigs.fa")
  tuple sample, file("assembled/${sample}_contigs.fa") into contigs_prokka
  file("logs/shovill/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p shovill ALL_assembled logs/shovill
    log_file=logs/shovill/!{sample}.!{workflow.sessionId}.log
    err_file=logs/shovill/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    shovill --version >> $log_file

    shovill --cpu !{task.cpus} \
      --ram !{task.memory} \
      --outdir shovill/!{sample} \
      --R1 !{reads[0]} \
      --R2 {reads[1]} \
      2>> $err_file | \
      tee $log_file
    cp shovill/!{sample}/contigs.fa assembled/!{sample}_contigs.fa

    exit
  '''
}

process mash_sketch {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  when:
  params.mash

  input:
  set val(sample), file(reads) from clean_reads_mash

  output:
  tuple sample, file("mash/${sample}.msh") into mash_sketch_files
  file("logs/mash_sketch/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p mash logs/mash_sketch
    log_file=logs/mash_sketch/!{sample}.!{workflow.sessionId}.log
    err_file=logs/mash_sketch/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file

    cat !{reads} | mash sketch -m 2 -o mash/!{sample} - 2>> $err_file | tee $log_file


    exit
  '''
}

process mash_dist {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus params.maxcpus

  when:
  params.mash

  input:
  set val(sample), file(msh) from mash_sketch_files

  output:
  tuple sample, file("mash/${sample}_mashdist.txt")
  tuple sample, env(genus), env(species), env(genome_size) into mash_results, mash_prokka
  file("logs/mash_dist/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p mash logs/mash_dist
    log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
    err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file

    mash dist -p !{task.cpus} -v 0 /db/RefSeqSketchesDefaults.msh !{msh} | sort -gk3 > {output.file} 2>> $err_file

    species="NA"

    grep $genus $species !{genome} | cut -f 2
    exit
  '''
}

contigs_prokka
  .join(mash_prokka, remainder: true, by:0)
  .set { for_prokka }

process prokka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  when:
  params.prokka

  input:
  set val(sample), file(contigs), val(genus), val(species) from for_prokka

  output:
  file("prokka/${sample}/${sample}.gff")
  file("gff/${sample}.gff")
  file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p mash logs/mash_dist
    log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
    err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file

    prokka --cpu !{task.cpus} \
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

    exit
  '''
}

process prokka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  when:
  params.quast

  input:
  set val(sample), file(contigs) from contigs_quast

  output:
  file("quast/${sample}/report.tsv")
  file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

  shell:
    '''
    mkdir -p mash logs/mash_dist
    log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
    err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file

    quast.py --version >> $log_file
    quast.py {input} --output-dir quast/{wildcards.sample} --threads {threads} 2>> $err_file | tee -a $log_file
  '''
}

process CG_pipeline_shuffle {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  when:
  params.cg_pipeline

  input:
  set val(sample), file(reads) from clean_reads_cg

  output:
  file("Sequencing_reads/shuffled/${sample}_raw_shuffled.fastq.gz")
  file("logs/cg_shuffle/${sample}.${workflow.sessionId}.{log,err}")

  shell:
    '''
    mkdir -p cg_shuffle logs/cg_shuffle
    log_file=logs/cg_shuffle/!{sample}.!{workflow.sessionId}.log
    err_file=logs/cg_shuffle/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    run_assembly_shuffleReads.pl -gz !{reads}  2>> $err_file > cg_shuffle/!{sample}_raw_shuffled.fastq.gz
    '''
}

process CG_pipeline {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus params.maxcpus

  when:
  params.cg_pipeline

  input:
  file("cg-pipeline/${sample}.gc_pipeline.txt")
  output:
  file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

  shell:
    '''
      run_assembly_readMetrics.pl {input.shuffled_fastq} --fast --numcpus {threads} -e $genome_length 2>> {output.log}.err
    '''
}

process seqsero2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  when:
  params.seqsero2 && salmonella_flag

  input:

  output:
  tuple sample, env(species) into seqsero2_results
  file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

  shell:
    '''
    '''
}

rule seqsero:
    input:
        rules.seqyclean.output.read1,
        rules.seqyclean.output.read2
    output:
        file="SeqSero/{sample}/Seqsero_result.txt",
        final="SeqSero/{sample}.Seqsero_result.txt",
        log=temp("logs/seqsero/{sample}")
    threads:
        1
    singularity:
        "docker://staphb/seqsero:1.0.1"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        rm -R SeqSero/{wildcards.sample} || true
        SeqSero.py -m 2 -d SeqSero/{wildcards.sample} -i {input} #2>> {output.log}.err >> {output.log}.log || true
        cp {output.file} {output.final} || true
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
            '''
              mkdir -p mash logs/mash_dist
              log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
              err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

              # time stamp + capturing tool versions
              date | tee -a $log_file $err_file > /dev/null
              echo "mash version: $(mash --version | head -n 1 )" >> $log_file
            '''
        }
rule abricate:
    input:
        rules.shovill.output.file
    output:
        file="abricate_results/{database}/{database}.{sample}.out.tab",
        log=temp("logs/abricate/{sample}.{database}")
    threads:
        5
    singularity:
        "docker://staphb/abricate:0.8.13s"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        abricate --version >> {output.log}.log
        abricate --list >> {output.log}.log
        abricate --db {wildcards.database} --threads {threads} --minid 80 --mincov 80 {input} > {output.file} 2>> {output.log}.err || true
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
            '''
              mkdir -p mash logs/mash_dist
              log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
              err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

              # time stamp + capturing tool versions
              date | tee -a $log_file $err_file > /dev/null
              echo "mash version: $(mash --version | head -n 1 )" >> $log_file
            '''
        }
rule abricate_summary:
    input:
        expand("abricate_results/{database}/{database}.{sample}.out.tab", sample=SAMPLE, database=DATABASE),
    output:
        file="abricate_results/summary/{database}.abricate_summary.txt",
        log=temp("logs/abricate/{database}_summary")
    threads:
        1
    singularity:
        "docker://staphb/abricate:0.8.13s"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        abricate --version >> {output.log}.log
        abricate --summary abricate_results*/{wildcards.database}/{wildcards.database}*tab > {output.file} 2>> {output.log}.err || true
        cat abricate_results/{wildcards.database}/*out.tab > abricate_results/{wildcards.database}.all.txt
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }

rule bwa_index:
    input:
        rules.shovill.output.file
    output:
        index="shovill_result/{sample}/contigs.fa.sa",
        log=temp("logs/bwa/{sample}_index")
    singularity:
        "docker://staphb/shovill:1.0.4"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "bwa $(bwa 2>&1 | grep Version )" >> {output.log}.log
        bwa index {input} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
            '''
              mkdir -p mash logs/mash_dist
              log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
              err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

              # time stamp + capturing tool versions
              date | tee -a $log_file $err_file > /dev/null
              echo "mash version: $(mash --version | head -n 1 )" >> $log_file
            '''
        }
rule blastn:
    input:
        rules.shovill.output.file
    output:
        tsv="blast/{sample}.tsv",
        log=temp("logs/blastn/{sample}")
    threads:
        10
    singularity:
        "docker://ncbi/blast:2.9.0"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        blastn -version >> {output.log}.log
        echo "The blastdb location is $BLASTDB" >> {output.log}.log
        blastn -query {input} \
            -out {output.tsv} \
            -num_threads {threads} \
            -db /blast/blastdb/nt \
            -outfmt '6 qseqid staxids bitscore std' \
            -max_target_seqs 10 \
            -max_hsps 1 \
            -evalue 1e-25 2>> {output.log}.log | tee -a {output.log}.err || true
        touch {output}
        """

        process bwa {
          publishDir "${params.outdir}", mode: 'copy', pattern: "logs/bwa/*.{log,err}"
          tag "${sample}"
          echo false
          cpus params.maxcpus

          when:
          params.aligner == 'bwa'

          input:
          set val(sample), file(reads), file(reference_genome) from clean_reads_bwa

          output:
          tuple sample, file("aligned/${sample}.sam") into bwa_sams
          file("logs/bwa/${sample}.${workflow.sessionId}.{log,err}")
          tuple sample, env(bwa_version) into bwa_version

          shell:
          '''
            mkdir -p aligned logs/bwa
            log_file=logs/bwa/!{sample}.!{workflow.sessionId}.log
            err_file=logs/bwa/!{sample}.!{workflow.sessionId}.err

            # time stamp + capturing tool versions
            date | tee -a $log_file $err_file > /dev/null
            echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
            bwa_version="bwa : "$(bwa 2>&1 | grep Version)

            # index the reference fasta file
            bwa index !{reference_genome}

            # bwa mem command
            bwa mem -t !{task.cpus} !{reference_genome} !{reads} 2>> $err_file > aligned/!{sample}.sam
          '''
        }

rule bwa:
    input:
        contig=rules.shovill.output.file,
        read1=rules.seqyclean.output.read1,
        read2=rules.seqyclean.output.read2,
        index=rules.bwa_index.output.index,
        test=rules.blastn.output.tsv
    threads:
        48
    output:
        bam="bwa/{sample}.sorted.bam",
        bai="bwa/{sample}.sorted.bam.bai",
        log=temp("logs/bwa/{sample}")
    singularity:
        "docker://staphb/shovill:1.0.4"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "bwa $(bwa 2>&1 | grep Version )" >> {output.log}.log
        samtools --version >> {output.log}.log
        if [ -s "{input.test}" ]
        then
            bwa mem -t {threads} {input.contig} {input.read1} {input.read2} 2>> {output.log}.err | \
            samtools sort -o {output.bam} 2>> {output.log}.err > {output.bam} || true
        else
            echo "BLASTDB not found" | tee -a {output.log}.log {output.log}.err
        fi
        samtools index {output.bam} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """
        process sort {
          publishDir "${params.outdir}", mode: 'copy'
          tag "${sample}"
          echo false
          cpus params.maxcpus

          input:
          set val(sample), file(sam) from sams

          output:
          tuple sample, file("aligned/${sample}.sorted.bam") into pre_trim_bams, pre_trim_bams2
          tuple sample, file("aligned/${sample}.sorted.bam"), file("aligned/${sample}.sorted.bam.bai") into pre_trim_bams_bamsnap
          file("logs/sort/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
            mkdir -p aligned logs/sort
            log_file=logs/sort/!{sample}.!{workflow.sessionId}.log
            err_file=logs/sort/!{sample}.!{workflow.sessionId}.err

            # time stamp + capturing tool versions
            date | tee -a $log_file $err_file > /dev/null
            samtools --version >> $log_file

            samtools sort !{sam} 2>> $err_file | \
              samtools view -F 4 -o aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file

            # indexing the bams
            samtools index aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file
          '''
        }

process prokka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  echo false
  cpus 1

  when:
  params.prokka

  input:

  output:
  file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

  shell:
    '''
      mkdir -p mash logs/mash_dist
      log_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.log
      err_file=logs/mash_dist/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      echo "mash version: $(mash --version | head -n 1 )" >> $log_file
    '''
}

rule blobtools_create:
    input:
        contig=rules.shovill.output.file,
        blast=rules.blastn.output.tsv,
        bam=rules.bwa.output.bam
    output:
        cov="blobtools/{sample}.{sample}.sorted.bam.cov",
        json="blobtools/{sample}.blobDB.json",
        log=temp("logs/blobtools/{sample}_create")
    threads:
        1
    singularity:
        "docker://chrishah/blobtools:v1.1.1"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "blobtools version $(blobtools -v)" >> {output.log}.log
        blobtools create -o blobtools/{wildcards.sample} -i {input.contig} -b {input.bam} -t {input.blast} 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }

rule blobtools_view:
    input:
        rules.blobtools_create.output.json,
    output:
        txt="blobtools/{sample}.blobDB.table.txt",
        log=temp("logs/blobtools/{sample}_view")
    threads:
        1
    singularity:
        "docker://chrishah/blobtools:v1.1.1"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "blobtools version $(blobtools -v)" >> {output.log}.log
        blobtools view -i {input} -o blobtools/ 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }

rule blobtools_plot:
    input:
        table=rules.blobtools_view.output.txt,
        json=rules.blobtools_create.output.json
    output:
        png="blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.bam0.png",
        covpng="blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.read_cov.bam0.png",
        txt="blobtools/{sample}.blobDB.json.bestsum.species.p8.span.100.blobplot.stats.txt",
        log=temp("logs/blobtools/{sample}_plot")
    threads:
        1
    singularity:
        "docker://chrishah/blobtools:v1.1.1"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        echo "blobtools version $(blobtools -v)" >> {output.log}.log
        blobtools plot -i {input.json} -o blobtools/ -r species --format png 2>> {output.log}.err | tee -a {output.log}.log || true
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }

rule mlst:
    input:
        expand("ALL_assembled/{sample}_contigs.fa", sample=SAMPLE)
    output:
        file="mlst/mlst.txt",
        log=temp("logs/mlst/mlst")
    threads:
        1
    singularity:
        "docker://staphb/mlst:2.17.6-cv1"
    shell:
        """
        date | tee -a {output.log}.log {output.log}.err
        mlst --version >> {output.log}.log
        mlst ALL_assembled/*_contigs.fa > mlst/mlst.txt 2>> {output.log}.err || true
        touch {output}
        """
        process prokka {
          publishDir "${params.outdir}", mode: 'copy'
          tag "$sample"
          echo false
          cpus 1

          when:
          params.prokka

          input:

          output:
          file("logs/prokka/${sample}.${workflow.sessionId}.{log,err}")

          shell:
          '''
          '''
        }

process bbduk {
  tag           "${sample}"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/bbtools:38.98'
  maxForks      10

  input:
  tuple val(sample), file(reads)

  output:
  tuple val(sample), file("bbduk/${sample}_rmphix_R{1,2}.fastq.gz"),  emit: fastq
  path "bbduk/*",                                                     emit: files
  path "bbduk/${sample}.phix.stats.txt",                              emit: stats
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log",    emit: log
  tuple val(sample), env(phix_reads),                                 emit: phix_reads

  shell:
  '''
    mkdir -p bbduk logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    bbduk.sh --version >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    bbduk.sh !{params.bbduk_options} \
      in1=!{reads[0]} \
      in2=!{reads[1]} \
      out1=bbduk/!{sample}_rmphix_R1.fastq.gz \
      out2=bbduk/!{sample}_rmphix_R2.fastq.gz \
      outm=bbduk/!{sample}.matched_phix.fq \
      ref=/opt/bbmap/resources/phix174_ill.ref.fa.gz \
      stats=bbduk/!{sample}.phix.stats.txt \
      threads=!{task.cpus} \
      | tee -a $log_file

    phix_reads=$(grep Matched bbduk/!{sample}.phix.stats.txt | cut -f 2)
  '''
}

process bbmap {
  tag           "${sample}"
  label         "maxcpus"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/bbtools:38.98'
  maxForks      10

  input:
  tuple val(sample), file(fastq), file(contigs)

  output:
  tuple val(sample), file("bbmap/${sample}.mapped_sorted.bam*"),   emit: bam
  path "bbmap/*txt",                                               emit: stats
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    bbmap.sh --version >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    bbmap.sh !{params.bbmap_options} \
      in1=!{fastq[0]} \
      in2=!{fastq[1]} \
      out=bbmap/!{sample}.mapped.sam \
      ref=!{contigs} \
      covstats=bbmap/!{sample}.constats.txt \
      covhist=bbmap/!{sample}.covhist.txt \
      basecov=bbmap/!{sample}.basecov.txt \
      bincov=bbmap/!{sample}.bincov.txt \
      threads=!{task.cpus} \
      bamscript=bs.sh \
      | tee -a $log_file
      
    sh bs.sh | tee -a $log_file
  '''
}

process blobtools_create {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'chrishah/blobtools:v1.1.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA cpus 14
  //#UPHLICA memory 60.GB
  //#UPHLICA time '45m'
  
  input:
  tuple val(sample), file(contig), file(blastn), file(bam)

  output:
  tuple val(sample), file("blobtools/${sample}.blobDB.json")     , emit: json
  path "blobtools/${sample}.${sample}*.bam.cov"                  , emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools create !{params.blobtools_create_options} \
      -o blobtools/!{sample} \
      -i !{contig} \
      -b !{bam[0]} \
      -t !{blastn} \
      | tee -a $log_file
  '''
}

process blobtools_view {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'chrishah/blobtools:v1.1.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA cpus 3
  //#UPHLICA memory 1.GB
  //#UPHLICA time '10m'
  
  input:
  tuple val(sample), file(json)

  output:
  tuple val(sample), file("blobtools/${sample}.blobDB.table.txt"), emit: file
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools view !{params.blobtools_view_options} \
      -i !{json} \
      -o blobtools/ \
      | tee -a $log_file
  '''
}

process blobtools_plot {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'chrishah/blobtools:v1.1.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA cpus 3
  //#UPHLICA memory 1.GB
  //#UPHLICA time '10m'
  
  input:
  tuple val(sample), file(json)

  output:
  path "blobtools/${sample}.*"                                   , emit: files
  tuple val(sample), file("blobtools/${sample}_blobtools.txt")   , emit: results
  path "blobtools/${sample}_summary.txt"                         , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools plot !{params.blobtools_plot_options} \
      -i !{json} \
      -o blobtools/ \
      | tee -a $log_file

    grep "^# " blobtools/!{sample}*.stats.txt | \
      sed 's/# //g' | \
      awk '{print "sample\t" $0 }' > blobtools/!{sample}_summary.txt

    grep -v "^#" blobtools/!{sample}*.stats.txt | \
      sed 's/%//g' | \
      tr " " "_" | \
      awk -v sample=!{sample} '{if ($13 >= 5.0 ) print sample "\\t" $0}' | \
      tr " " "\\t" | \
      sort -k 14rn,14 >> blobtools/!{sample}_summary.txt

    grep -vw all blobtools/!{sample}_summary.txt > blobtools/!{sample}_blobtools.txt
  '''
}
process legsta {
  tag           "$meta.id"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'staphb/legsta:0.5.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'

  when:
  flag =~ 'found'

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(sample), file(contigs), val(flag)

  output:
  path "legsta/${sample}_legsta.csv"                             , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p legsta logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    legsta --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file
    
    cat !{contigs} | awk '{($1=$1); print $0}' > !{sample}.fasta.tmp
    mv !{sample}.fasta.tmp !{sample}.fasta

    legsta !{params.legsta_options} \
      !{sample}.fasta \
      --csv \
      | tee -a $log_file \
      | awk -v sample=!{sample} '{print sample "," $0 }' \
      | sed '0,/!{sample}/{s/!{sample}/sample/}' > legsta/!{sample}_legsta.csv
  '''
}
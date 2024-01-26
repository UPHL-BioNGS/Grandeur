process kleborate {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/kleborate:2.4.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(contig), val(flag), file(script)

  output:
  path "kleborate/${sample}_results.tsv"                         , emit: collect, optional: true
  path "kleborate/${sample}_results.txt"                         , emit: result
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p kleborate logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    kleborate --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    kleborate !{params.kleborate_options} \
      -o kleborate/!{sample}_results.txt \
      -a !{contig} \
      | tee -a $log_file

    python3 !{script} kleborate/!{sample}_results.txt kleborate/!{sample}_results.tsv kleborate !{sample}
  '''
}

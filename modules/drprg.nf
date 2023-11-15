process drprg {
  tag           "${sample}"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/drprg:0.1.1--h5076881_1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'

  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(contigs), val(flag)

  output:
  path "drprg/${sample}_drprg.csv"                             , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p drprg logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    drprg --version >> $log_file
    drprg index --list >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file
    
    drprg predict !{params.drprg_options} \
      !{contigs} \
      -x /drprg/mtb/mtb \
      -i !{contigs} \
      -o drprg/ \
      | tee =a $log_file

    exit 1
  '''
}
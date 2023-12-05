process phytreeviz {
  tag           "${analysis}"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/phytreeviz:0.1.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA cpus 14
  //#UPHLICA memory 60.GB
  //#UPHLICA time '24h'
  
  input:
  tuple val(analysis), file(newick)

  output:
  path "phytreeviz/${analysis}_tree.png",                                            emit: for_multiqc
  path "logs/${task.process}/${analysis}.${task.process}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p phytreeviz logs/!{task.process}
    log_file=logs/!{task.process}/!{analysis}.!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    phytreeviz -v >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    phytreeviz !{params.phytreeviz_options} \
        -i !{newick} \
        -o phytreeviz/!{analysis}_tree.png
  '''
}

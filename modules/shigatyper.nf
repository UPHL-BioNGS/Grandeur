process shigatyper {
  tag           "${sample}"
  label         "medcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/shigatyper:2.0.5'
  stageInMode   'copy'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  //#UPHLICA memory 26.GB
  //#UPHLICA cpus 7
  //#UPHLICA time '10m'
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(input), val(flag), file(script)

  output:
  path "shigatyper/${sample}_shigatyper.tsv",      optional: true, emit: files
  path "shigatyper/${sample}_shigatyper-hits.tsv", optional: true, emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p shigatyper logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    shigatyper --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    shigatyper !{params.shigatyper_options} \
      --SE !{input} \
      --name !{sample} \
      | tee -a $log_file

    python3 !{script} !{sample}-hits.tsv shigatyper/!{sample}_shigatyper-hits.tsv shigatyper !{sample}

    if [ -f "!{sample}.tsv" ] ; then cp !{sample}.tsv shigatyper/!{sample}_shigatyper.tsv ; fi
  '''
}

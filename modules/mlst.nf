process mlst {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mlst:2.23.0-2023-07'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  tuple val(sample), file(contig), file(script)

  output:
  path "mlst/${sample}_mlst.tsv"                                 , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p mlst logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    mlst --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    mlst !{params.mlst_options} \
      --threads !{task.cpus} \
      !{contig} \
      > !{sample}_mlst.txt

    python3 !{script} !{sample}_mlst.txt mlst/!{sample}_mlst.tsv mlst !{sample}
  '''
}

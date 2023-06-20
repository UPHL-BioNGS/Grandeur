process snp_matrix {
  tag           "Visualizing SNP matrix"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/uphl/seaborn:0.12.2'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'
  
  input:
  file(matrix)
  file(script)

  output:
  path "snp-dists/SNPmatrix.pdf"                                       , emit: snp_matrix
  path "snp-dists/SNPmatrix_mqc.png"                                   , emit: for_multiqc
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p snp-dists logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    python !{script}

    mv SNPmatrix.pdf snp-dists/.
    mv SNPmatrix_mqc.png snp-dists/.
  '''
}

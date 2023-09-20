process srst2 {
  tag           "${sample}"
  label         "medcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/srst2:0.2.0-vcholerae'
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
  tuple val(sample), file(contigs), val(flag), file(script)

  output:
  path "srst2/${sample}*results.txt",      optional: true, emit: files
  path "srst2/*"
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p srst2 logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    srst2 --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    srst2 !{params.srst2_options} \
        --input_se !{contigs} \
        --gene_db /vibrio-cholerae-db/vibrio_230224.fasta \
        --read_type f \
        --output srst2/!{sample} \
        --threads !{task.cpus}

    #python3 !{script} input output srst2 !{sample}

    #if [ -f output ] ; then cp output srst2/!{sample}_srst2.csv ; fi
  '''
}

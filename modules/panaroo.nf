process panaroo {
  tag           "Core Genome Alignment"
  label         'maxcpus'
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/panaroo:1.3.4--pyhdfd78af_0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'hicpu-small'
  //#UPHLICA cpus 15
  //#UPHLICA memory 30.GB
  //#UPHLICA time '10m'
  
  input:
  file(contigs)

  output:
  path "panaroo/*"                                                                         , emit: files
  tuple path("panaroo/core_gene_alignment.aln"), path("panaroo/gene_presence_absence.Rtab"), emit: core_gene_alignment
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"                    , emit: log_files

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    panaroo --version >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    panaroo !{params.panaroo_options} \
        -t !{task.cpus} \
        -o panaroo \
        -i !{contigs} \
        -a core \
        | tee -a $log_file
  '''
}

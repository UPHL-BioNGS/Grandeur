process panaroo {
  tag           "Core Genome Alignment"
  label         "process_high"
  label         'maxcpus'
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/panaroo:1.3.4'
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
path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: '--clean-mode strict --remove-invalid-genes'
    def prefix = task.ext.prefix ?: "${meta.id}"
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

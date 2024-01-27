process quast {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/quast:5.2.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-xlarge'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 14
  //#UPHLICA time '10m'
  
  input:
  tuple val(meta), file(contigs)

  output:
  path "quast/*"                                                     , emit: files
  path "quast/*_quast_report.tsv"                    , optional: true, emit: for_multiqc
  tuple val(meta), file("quast/*_quast_report.tsv"), optional: true, emit: results
  path "quast/*/transposed_report.tsv"               , optional: true, emit: collect
  path "logs/${task.process}/*.log"            , emit: log
  path "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p ${task.process} logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    quast.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    quast.py ${params.quast_options} \
      ${contigs} \
      --output-dir quast/${prefix} \
      --threads ${task.cpus} \
      | tee -a $log_file

    if [ -f "quast/${prefix}/report.tsv" ] ; then cp quast/${prefix}/report.tsv quast/${prefix}_quast_report.tsv ; fi

    if [ -f "quast/${prefix}/transposed_report.tsv" ]
    then
      head -n 1 quast/${prefix}/transposed_report.tsv | awk '{print "prefix\\t" $0 }' > quast/${prefix}/transposed_report.tsv.tmp
      tail -n 1 quast/${prefix}/transposed_report.tsv | awk -v prefix=${prefix} '{print prefix "\\t" $0}' >> quast/${prefix}/transposed_report.tsv.tmp
      mv quast/${prefix}/transposed_report.tsv.tmp quast/${prefix}/transposed_report.tsv
    fi
  """
}

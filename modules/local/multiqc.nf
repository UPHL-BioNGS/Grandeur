process multiqc {
  tag           "multiqc"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(input)
  file(script)

  output:
  path "multiqc/multiqc_report.html", optional: true, emit: report
  path "multiqc/multiqc_data/*"     , optional: true, emit: data_folder
  path "logs/${task.process}/*.log" , emit: log_files

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: ''
  """
    mkdir -p multiqc quast logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    for quast_file in \$(ls *_quast_report.tsv)
    do
      sample=\$(echo \$quast_file | sed 's/_quast_report.tsv//g' | head -n 1 )
      mkdir -p quast/\$sample
      mv \$quast_file quast/\$sample/report.tsv
    done

    python ${script}

    multiqc ${args} \
      --outdir multiqc \
      --cl_config "prokka_fn_snames: True"  \
      . \
      | tee -a \$log_file
  """
}

process versions {
  tag           "version dump"
  label         "process_single"
  // no publishDir
  container     'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(input)
  file(versions_script)

  output:
  path "multiqc/multiqc_report.html", optional: true, emit: report
  path "multiqc/multiqc_data/*"     , optional: true, emit: data_folder
  path "logs/${task.process}/*.log" , emit: log_files
  path  "versions.yml"              , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  """

    exit 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS

    python ${versions_script}

    exit 1
  """
}
process multiqc {
  tag           "multiqc"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/multiqc:1.19'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(input)

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

    multiqc ${args} \
      --outdir multiqc \
      --cl-config "prokka_fn_snames: True"  \
      . \
      | tee -a \$log_file
  """
}

process versions {
  tag           "extracting versions"
  label         "process_single"
  // no publishDir
  container     'staphb/multiqc:1.19'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(input)
  file(versions_script)

  output:
  path "software_versions_mqc.yml", emit: for_multiqc
  path "software_versions.yml", emit: yml

  when:
  task.ext.when == null || task.ext.when

  shell:
  """
    cat <<-END_VERSIONS >> versions.yml
    "report:multiqc":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS

    python3 ${versions_script}
  """
}
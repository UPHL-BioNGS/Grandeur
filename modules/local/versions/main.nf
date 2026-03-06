process VERSIONS {
  tag           "extracting versions"
  label         "process_single"
  container     'staphb/multiqc:1.33'

  input:
  file(input)
  file(versions_script)

  output:
  path "software_versions_mqc.yml", emit: for_multiqc
  path "software_versions.yml", emit: yml

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    cat <<-END_VERSIONS >> versions.yml
    "REPORT:MULTIQC":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS

    python3 ${versions_script}
  """
}

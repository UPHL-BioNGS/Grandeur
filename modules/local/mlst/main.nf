process MLST {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/mlst:2.35.0'

  input:
  tuple val(meta), file(contig)

  output:
  tuple val(meta), file("mlst/*_mlst.txt"), emit: files, optional: true
  path "*_mlst_summary.txt", emit: collect, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p mlst 

    # there are file permission issues sometimes
    cat ${contig} > input_${prefix}.fasta

    mlst ${args} \
      --full \
      --threads ${task.cpus} \
      input_${prefix}.fasta | \
      tr ' ' '_' \
      > mlst/${prefix}_mlst.txt

    if [ -f "mlst/${prefix}_mlst.txt" ]
    then
      head -n 1  mlst/${prefix}_mlst.txt | awk '{print "sample\\t" \$0 }' > ${prefix}_mlst_summary.txt
      tail -n 1  mlst/${prefix}_mlst.txt| awk -v sample=${prefix} '{print sample "\\t" \$0}' >> ${prefix}_mlst_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
        scheme_db_date: 2026-01-13
    END_VERSIONS
  """
}

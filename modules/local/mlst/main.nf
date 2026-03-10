process MLST {
  tag           "${meta.id}"
  label         "process_medium"
  container     'staphb/mlst:2.32.2'

  input:
  tuple val(meta), file(contig)

  output:
  tuple val(meta), file("mlst/*_mlst.txt"), emit: files, optional: true
  path "*_mlst_summary.txt", emit: collect
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
      --threads ${task.cpus} \
      input_${prefix}.fasta | \
      tr ' ' '_' \
      > mlst/${prefix}_mlst.txt

    echo -e "sample\\tfilename\\tmatching PubMLST scheme\\tST\\tID1\\tID2\\tID3\\tID4\\tID5\\tID6\\tID7\\tID8\\tID9\\tID10\\tID11\\tID12\\tID13\\tID14\\tID15" > ${prefix}_mlst_summary.txt
    cat mlst/${prefix}_mlst.txt | awk '{print "${prefix}\\t" \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11 "\\t" \$12 "\\t" \$13 "\\t" \$14 "\\t" \$15 "\\t" \$16 "\\t" \$17 "\\t" \$18 }' >> ${prefix}_mlst_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
        scheme_db_date: 2026-01-13
    END_VERSIONS
  """
}

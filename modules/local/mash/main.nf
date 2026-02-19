process MASH {
  tag        "${meta.id}"
  label      "process_medium"
  container  'staphb/mash:2.3'

  input:
  tuple val(meta), file(reads), file(reference)

  output:
  path "mash/*.mashdist.txt",                       emit: mashdist
  tuple val(meta), file("mash/*.summary.mash.csv"), emit: results
  path "*err",                                      emit: mash_err
  path "logs/${task.process}/*.log",                emit: log
  path "versions.yml",                              emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
    def args        = task.ext.args        ?: "-v 0 -d 0.25"
    def args_sketch = task.ext.sketch_args ?: "-s 1000 -k 21" 
    def is_fastq    = (reads instanceof List) || reads.name.toString().matches('.*\\.(fastq|fq)(\\.gz)?$')
    def sketch_opts = is_fastq ? "-m 2 ${args_sketch}" : "-m 1 ${args_sketch}"    
    def better_ref  = reference.toString().contains("input") ? "/db/RefSeqSketchesDefaults.msh" : reference
    def org_cmd     = reference.toString().contains("input") ? "echo \$line | cut -f 8 -d - | cut -f 1,2 -d _ | cut -f 1 -d ." : "echo \$line | cut -f 1,2 -d _ | cut -f 1 -d ."
    def prefix      = task.ext.prefix      ?: "${meta.id}"
    def input_files = (reads instanceof List) ? reads.join(" ") : reads

    """
    mkdir -p mash logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # mash sketch, err_file contains estimated genome size
    cat ${input_files} | \
      mash sketch ${sketch_opts} \
      -p ${task.cpus} \
      -o ${prefix} - \
      2>> mash_${prefix}.err | tee -a \$log_file

    # mash dist
    mash dist ${args} \
      -p ${task.cpus} \
      ${better_ref} \
      ${prefix}.msh | \
      sort -gk3 \
      > mash/${prefix}.mashdist.txt

    echo "sample,reference,query,mash-distance,P-value,matching-hashes,organism" > mash/${prefix}.summary.mash.csv

    while read line
    do
      organism=\$(${org_cmd})
      echo \$line | \
        sed 's/,//g' | \
        awk -v sample=${prefix} \
        -v org=\$organism \
        '{print sample "," \$1 "," \$2 "," \$3 "," \$4 "," \$5 "," org}' >> mash/${prefix}.summary.mash.csv
    done < mash/${prefix}.mashdist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mash: \$( mash --version )
    END_VERSIONS
  """
}
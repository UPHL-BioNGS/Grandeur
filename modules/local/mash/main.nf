process MASH {
  tag        "${meta.id}"
  label      "process_medium"
  container  'staphb/mash:2.3'

  input:
  tuple val(meta), file(reads), file(reference)

  output:
  tuple val(meta), file("mash/*.mashdist.txt"),   emit: mashdist, optional: true
  tuple val(meta), file("mash/*.mashscreen.txt"), emit: screen, optional: true
  path "mash/*_mashdist_summary.csv",             emit: results, optional: true
  path "mash/*_mashscreen_summary.csv",           emit: screen_results, optional: true
  path "*err",                                    emit: mash_err, optional: true
  path "logs/${task.process}/*.log",              emit: log
  path "versions.yml",                            emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
    def args        = task.ext.args        ?: "-v 0 -d 0.15"
    def screen_args = task.ext.screen_args ?: "-v 0 -w"
    def args_sketch = task.ext.sketch_args ?: "-s 1000 -k 21" 
    def is_fastq    = (reads instanceof List) || reads.name.toString().matches('.*\\.(fastq|fq)(\\.gz)?$')
    def sketch_opts = is_fastq ? "-m 2 ${args_sketch}" : "-m 1 ${args_sketch}"    
    def better_ref  = reference.toString().contains("input") ? "/db/RefSeqSketchesDefaults.msh" : reference
    def org_cmd     = reference.toString().contains("input") ? "echo \$line | cut -f 8 -d - | cut -f 1,2 -d _ | cut -f 1 -d ." : "echo \$line | cut -f 1,2 -d _ | cut -f 1 -d ."
    def org_cmd2    = reference.toString().contains("input") ? "echo \$line | awk '{print \$5}' | cut -f 8 -d - | cut -f 1,2 -d _ | cut -f 1 -d ." : "echo \$line | awk '{print \$5}' | cut -f 1,2 -d _ | cut -f 1 -d ."
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

    echo "sample,reference,query,mash-distance,P-value,matching-hashes,organism" > mash/${prefix}_mashdist_summary.csv

    while read line
    do
      organism=\$(${org_cmd})
      echo \$line | \
        sed 's/,//g' | \
        awk -v sample=${prefix} \
        -v org=\$organism \
        '{print sample "," \$1 "," \$2 "," \$3 "," \$4 "," \$5 "," org}' >> mash/${prefix}_mashdist_summary.csv
    done < mash/${prefix}.mashdist.txt

    # mash screen
    mash screen ${screen_args} \
      -p ${task.cpus} \
      ${better_ref} \
      ${reads[0]} \
      > mash/${prefix}.mashscreen.txt

    echo "sample,identity,shared-hashes,median-multiplicity,p-value,query-ID,organism" > mash/${prefix}_mashscreen_summary.csv

    while read line
    do
      organism2=\$(${org_cmd2})
      echo \$line | \
        sed 's/,//g' | \
        awk -v sample=${prefix} \
        -v org=\$organism2 \
        '{print sample "," \$1 "," \$2 "," \$3 "," \$4 "," \$5 "," org}' >> mash/${prefix}_mashscreen_summary.csv
    done < mash/${prefix}.mashscreen.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mash: \$( mash --version )
    END_VERSIONS
  """
}
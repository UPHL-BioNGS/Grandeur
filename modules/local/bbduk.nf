process bbduk {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/bbtools:39.01'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("bbduk/*_rmphix_R{1,2}.fastq.gz"),  emit: fastq
  path "bbduk/*",                                           emit: files
  path "bbduk/*.phix.stats.txt",                            emit: stats
  path "logs/${task.process}/*.log",  emit: log
  path "versions.yml",                                      emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: 'k=31 hdist=1'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p bbduk logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    bbduk.sh ${args} \
      in1=${reads[0]} \
      in2=${reads[1]} \
      out1=bbduk/${prefix}_rmphix_R1.fastq.gz \
      out2=bbduk/${prefix}_rmphix_R2.fastq.gz \
      outm=bbduk/${prefix}.matched_phix.fq \
      ref=/opt/bbmap/resources/phix174_ill.ref.fa.gz \
      stats=bbduk/${prefix}.phix.stats.txt \
      threads=${task.cpus} \
      | tee -a \$log_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      bbduk: "\$(bbduk.sh --version 2>&1 | grep -v java | grep version | awk '{print \$NF}')"
    END_VERSIONS
  """
}

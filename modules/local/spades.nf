process spades {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    path: params.outdir, mode: 'copy', pattern: 'logs/*/*log'
  publishDir    path: params.outdir, mode: 'copy', pattern: 'spades/*'
  publishDir    path: params.outdir, mode: 'copy', pattern: 'spades/*/*'
  publishDir    path: params.outdir, mode: 'copy', pattern: 'contigs/*'
  container     'staphb/spades:3.15.5'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '5h'
  
  input:
  tuple val(meta), file(reads)

  output:
  path "spades/*/*", emit: files
  tuple val(meta), file("contigs/*_contigs.fa"), optional: true,  emit: contigs
  tuple val(meta), file("contigs/*_contigs.fa"), file(reads), optional: true,  emit: reads_contigs
  path "logs/${task.process}/*.log", emit: log
  path  "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '--isolate'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p spades contigs logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    spades.py ${args} \
      -1 ${reads[0]} \
      -2 ${reads[1]} \
      --threads ${task.cpus} \
      -o spades/${prefix} \
      | tee -a \$log_file

    if [ -f "spades/${prefix}/contigs.fasta" ] ; then cp spades/${prefix}/contigs.fasta contigs/${prefix}_contigs.fa ; fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
  """
}

process prokka {
  tag           "$meta.id"
  label         "process_high"
  label         "maxcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/prokka:1.14.6'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'himem-small'
  //#UPHLICA memory 60.GB
  //#UPHLICA cpus 7
  //#UPHLICA time '24h'

  input:
  tuple val(meta), file(contigs), val(organism)

  output:
  path "prokka/*/*", emit: prokka_files
  path "prokka/*/*.txt", emit: for_multiqc
  path "gff/*.gff", emit: gffs, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: '--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB'
  def prefix = task.ext.prefix ?: "${meta.id}"
  if (organism $= null) {
    """
      mkdir -p prokka gff logs/${task.process}
      log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

      # time stamp + capturing tool versions
      date > $log_file
      echo "container : ${task.container}" >> $log_file
      prokka -v >> $log_file
      echo "Nextflow command : " >> $log_file
      cat .command.sh >> $log_file

      prokka ${params.prokka_options} \
        --cpu ${task.cpus} \
        --outdir prokka/${prefix} \
        --prefix ${prefix} \
        --genus ${organism[0]} \
        --species ${organism[1]} \
        --force ${contigs} \
        | tee -a $log_file

      cp prokka/${prefix}/${prefix}.gff gff/${prefix}.gff
    """
  } else {
    """
      mkdir -p prokka gff logs/${task.process}
      log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

      # time stamp + capturing tool versions
      date > $log_file
      echo "container : ${task.container}" >> $log_file
      prokka -v >> $log_file
      echo "Nextflow command : " >> $log_file
      cat .command.sh >> $log_file

      prokka ${params.prokka_options} \
        --cpu ${task.cpus} \
        --outdir prokka/${prefix} \
        --prefix ${prefix} \
        --force ${contigs} \
        | tee -a $log_file

      cp prokka/${prefix}/${prefix}.gff gff/${prefix}.gff
    """
  }
}

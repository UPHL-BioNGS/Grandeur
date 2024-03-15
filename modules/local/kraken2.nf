process kraken2 {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/kraken2:2.1.3'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'
  
  input:
  tuple val(meta), file(fastq), path(kraken2_db)

  output:
  path "kraken2/*_kraken2_report.txt",                          emit: for_multiqc
  path "kraken2/*",                                             emit: files
  path "logs/${task.process}/*.log",                            emit: log
  tuple val(meta), file("kraken2/*_reads_summary_kraken2.csv"), emit: results
  path  "versions.yml",                                         emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p kraken2 logs/${task.process}
  log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

  kraken2 ${args} \
    --paired \
    --classified-out kraken2/${prefix}.cseqs#.fastq \
    --unclassified-out kraken2/${prefix}.useqs#.fastq \
    --output ${prefix}.kraken2.classifiedreads.txt \
    --threads ${task.cpus} \
    --db ${kraken2_db} \
    ${fastq[0]} ${fastq[1]} \
    --report kraken2/${prefix}_kraken2_report.txt \
    | tee -a \$log_file

  echo "Sample,Percentage of fragments,Number of fragments,Number of fragments assigned directly to this taxon,Rank code,NCBI taxonomic ID number,Scientific name" > kraken2/${prefix}_reads_summary_kraken2.csv
  cat kraken2/${prefix}_kraken2_report.txt | grep -w S | sed 's/,//g' | \
    awk -v sample=${prefix} '{ if (\$1 >= 5 ) print sample "," \$1 "," \$2 "," \$3 "," \$4 "," \$5 "," \$6 "_" \$7 }' | \
    sort >> kraken2/${prefix}_reads_summary_kraken2.csv

  fastq_check=\$(ls kraken2/* | grep fastq\$ | head -n 1)
  if [ -n "\$fastq_check" ] ; then gzip kraken2/*fastq ; fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
  END_VERSIONS
  """
}

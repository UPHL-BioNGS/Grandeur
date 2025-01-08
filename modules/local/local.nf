process CORE_GENOME_EVALUATION {
  tag         "Evaluating core genome"
  label       "process_single"
  container   'quay.io/biocontainers/pandas:1.5.2'

  input:
  tuple file(fasta), file(summary), file(script)

  output:
  tuple file(fasta), env("num_samples"), env("num_core_genes"), emit: evaluation
  path "core_genome_evaluation/core_genome_evaluation.csv", emit: for_multiqc
  path "logs/${task.process}/*.log"                       , emit: log_files

  script:
  """
    mkdir -p core_genome_evaluation logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    python ${script} | tee -a \$log_file

    num_samples=\$(wc -l core_genome_evaluation.csv | awk '{print \$1}' )
    num_core_genes=\$(cut -f 3 core_genome_evaluation.csv -d "," | tail -n 1 | cut -f 1 -d "." )
    cp core_genome_evaluation.csv core_genome_evaluation/core_genome_evaluation.csv
  """
}

process DOWNLOAD_SRA {
  tag           "${SRR}"
  label         "process_single"
  container     'quay.io/biocontainers/pandas:1.5.2'
  
  input:
  val(SRR)

  output:
  tuple val(SRR), file("reads/${SRR}_{1,2}.fastq.gz"), emit: fastq
  path "logs/${task.process}/*.log", emit: log

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    mkdir -p reads logs/${task.process}
    log_file=logs/${task.process}/${SRR}.${workflow.sessionId}.log

    echo "fasterq-dump failed. Attempting download from ENA" | tee -a \$log_file
      
    sra=${SRR}

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/\${sra:0:6}/0\${sra: -2}/${SRR}/${SRR}_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/\${sra:0:6}/0\${sra: -2 }/${SRR}/${SRR}_2.fastq.gz

    mv *fastq.gz reads/.
  """
}

process JSON_CONVERT {
  tag       "${meta.id}"
  label     "process_single"
  container 'quay.io/biocontainers/pandas:1.5.2'

  input:
  tuple val(meta), val(analysis), file(json), file(script)

  output:
  path "${analysis}/*_${analysis}*", emit: collect

  script:
  """
  mkdir -p ${analysis}

  python3 ${script} ${json} ${analysis}

  mv *${analysis}*tsv ${analysis}/.
  """
}

process MASH_ERR {
  tag           "${meta.id}"
  label         "process_single"
  container     'quay.io/biocontainers/pandas:1.5.2'
  
  input:
  tuple val(meta), file(error_file)

  output:
  path "mash_estimates.csv", emit: summary

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """

  genome_size=\$(grep "Estimated genome size:" ${error_file} | awk '{print \$NF}')
  coverage=\$(grep "Estimated coverage:"  ${error_file} | awk '{print \$NF}')

  echo "sample,mash_estimated_genome_size,mash_estimated_coverage" > mash_estimates.csv
  echo "${prefix},\$genome_size,\$coverage" >> mash_estimates.csv
  """
}

process MQC_PREP {
  tag           "prepping files"
  label         "process_single"
  container     'quay.io/biocontainers/pandas:1.5.2'
  
  input:
  file(input)
  file(script)

  output:
  path "*mqc*", emit: for_multiqc

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  python3 ${script}
  """
}

process NAMES {
  tag           "${meta.id}"
  label         "process_single"
  container     'quay.io/biocontainers/pandas:1.5.2'
  
  input:
  tuple val(meta), file(input)

  output:
  path "summary/*_names.csv", emit: collect

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def files  = input.join(" ")
  """
  mkdir -p summary

  echo "sample,file,version" > summary/${prefix}_names.csv
  echo "${prefix},${files},${workflow.manifest.version}" >> summary/${prefix}_names.csv
  """
}

process REFERENCES {
  tag       "Preparing references"
  label     "process_single"
  container 'quay.io/uphl/grandeur_ref:2024-06-26'

  output:
  path "ref/*", emit: fastas

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  mkdir ref

  cp /ref/* ref/.
  """
}

process SPECIES {
  tag           "Creating list of species"
  label         "process_single"
  container     'quay.io/biocontainers/pandas:1.5.2'
  
  input:
  file(results)

  output:
  path "datasets/species_list.txt", emit: species                                    

  when:
  task.ext.when == null || task.ext.when

  script:

  """
    mkdir -p datasets

    if [ -f "blobtools_species.txt" ]
    then 
      cut -f 2 blobtools_species.txt >> species.txt
    fi

    if [ -f "kraken2_summary.csv" ]
    then 
      cut -f 7 -d , kraken2_summary.csv >> species.txt
    fi

    if [ -f "mash_summary.csv" ]
    then
      cut -f 7 -d , mash_summary.csv | tail -n+2 >> species.txt
    fi

    grep -v no-hit species.txt | grep -v undef | grep -v name | grep "_" | sed 's/^_//g' | sort | uniq > datasets/species_list.txt
  """
}

process SUMMARY {
  tag           "Creating summary files"
  container     'quay.io/biocontainers/pandas:1.5.2'
  label         "process_single"
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(input)

  output:
  path "grandeur_summary.tsv"                 , emit: summary_tsv
  path "grandeur_summary.txt"                 , emit: summary_txt
  path "summary/grandeur_extended_summary.tsv", emit: extended_tsv
  path "summary/grandeur_extended_summary.txt", emit: extended_txt

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    mkdir -p summary

    python summary.py
  """
}

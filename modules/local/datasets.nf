process datasets_summary {
  tag           "${taxon}"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/uphl/datasets:16.3.0-2024-01-23'
  time          '1h'  
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
  
  input:
  tuple val(taxon), file(script)

  output:
  path "datasets/*_genomes.csv", emit: genomes, optional: true
  path "logs/${task.process}/*.${workflow.sessionId}.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: ''
  """
    mkdir -p datasets logs/${task.process}
    log_file=logs/${task.process}/${taxon}.${workflow.sessionId}.log 

    python3 ${script} ${taxon} ${params.datasets_max_genomes} \
    | tee -a \$log_file
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      "datasets: \$(datasets --version)"
    END_VERSIONS
  """
}

// It is faster if datasets can download the entire list at a time, but there is a 20 minute timeout for downloading.
// The '||' is to allow each genome to be downloaded on its own, which is longer overall but each genome should be less than 20 minutes.
process datasets_download {
  tag           "Downloading Genomes"
  // because there's no way to specify threads
  label         "process_medium"
  publishDir    path: "${params.outdir}", mode: 'copy', pattern: "logs/*/*log"
  container     'quay.io/uphl/datasets:16.3.0-2024-01-23'
  time          '5h'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  file(ids)

  output:
  path "genomes/*",    emit: genomes, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  """
    mkdir -p datasets genomes

    grep -h -v accession *csv | cut -f 1 -d , | sort | uniq > this_run.txt

    cat all_runs.txt this_run.txt | sort | uniq > id_list.txt

    ( datasets download genome accession --inputfile id_list.txt --filename ncbi_dataset.zip ; unzip -o ncbi_dataset.zip ) || \
    ( while read line ; do echo "Downloading \$line" ; datasets download genome accession \$line --filename dataset.zip ; unzip -o dataset.zip ; done < id_list.txt )

    fastas=\$(ls ncbi_dataset/data/*/*.fna )

    for fasta in \${fastas[@]}
    do
      echo "Copying \$fasta to genomes"
      accession=\$(echo \$fasta | cut -f 4 -d / | cut -f 1,2 -d _ )
      organism=\$(head -n 1 \$fasta | awk '{print \$2 "_" \$3 }' | sed 's/,//g' | sed 's/\\]//g' | sed 's/\\[//g' )
      cat \$fasta | sed 's/ /_/g' | sed 's/,//g' > genomes/\${organism}_\${accession}_ds.fna
      gzip genomes/\${organism}_\${accession}_ds.fna
    done  

    # removing MAGS
    rm -rf genomes/*:_*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      "datasets: \$(datasets --version)"
    END_VERSIONS
  """
}

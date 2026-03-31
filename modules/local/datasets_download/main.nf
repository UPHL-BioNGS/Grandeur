// It is faster if datasets can download the entire list at a time, but there is a 20 minute timeout for downloading.
// The '||' is to allow each genome to be downloaded on its own, which is longer overall but each genome should be less than 20 minutes.
process DATASETS_DOWNLOAD {
  tag           "Downloading Genomes"
  label         "process_medium"
  container     'staphb/ncbi-datasets:18.18.0'
  
  input:
  file(ids)

  output:
  path "genomes/*",    emit: genomes, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    mkdir -p datasets genomes

    grep -h -v accession *csv | cut -f 1 -d , | sort | uniq > this_run.txt

    cat id_list.txt this_run.txt | sort | uniq > accessions.txt

    ( datasets download genome accession --inputfile accessions.txt --filename ncbi_dataset.zip ; unzip -o ncbi_dataset.zip ) || \
    ( while read line ; do echo "Downloading \$line" ; datasets download genome accession \$line --filename dataset.zip ; unzip -o dataset.zip ; done < id_list.txt )

    fastas=\$(ls ncbi_dataset/data/*/*.fna )

    for fasta in \${fastas[@]}
    do
      echo "Copying \$fasta to genomes"
      accession=\$(echo \$fasta | cut -f 4 -d / | cut -f 1,2 -d _ )
      organism=\$(head -n 1 \$fasta | awk '{print \$2 "_" \$3 }' | sed 's/,//g' | sed 's/\\]//g' | sed 's/\\[//g' )
      cat \$fasta | sed 's/ /_/g' | sed 's/,//g' > genomes/\${organism}_\${accession}.fna
      gzip genomes/\${organism}_\${accession}.fna
    done  

    # removing MAGS
    rm -rf genomes/*:_*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      datasets: \$(datasets --version | awk '{print \$NF}')
    END_VERSIONS
  """
}

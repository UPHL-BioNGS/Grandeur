process DATASETS_SUMMARY {
  tag           "${taxon}"
  label         "process_single"
  container     'staphb/ncbi-datasets:16.35.0'

  input:
  tuple val(taxon), file(script)

  output:
  path "datasets/*_genomes.csv", emit: genomes, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args     = task.ext.args  ?: '--reference --mag exclude'
  def args2    = task.ext.args2 ?: '--annotated --assembly-level complete,scaffold --mag exclude'
  def fields   = task.ext.fields ?: 'accession,assminfo-refseq-category,assminfo-level,organism-name,assmstats-total-ungapped-len'
  def gen_spec = taxon.replaceAll(~/_.sp$/,"").split('_').join(" ")
  """
    mkdir -p datasets logs/!{task.process}

    datasets summary genome taxon "${gen_spec}" ${args}  --limit ${params.datasets_max_genomes} --as-json-lines | \
      dataformat tsv genome --fields ${fields} | \
      tee -a ${taxon}_genomes_rep.tsv

    datasets summary genome taxon "${gen_spec}" ${args2} --limit ${params.datasets_max_genomes} --as-json-lines | \
      dataformat tsv genome --elide-header --fields ${fields} | \
      tee -a ${taxon}_genomes_etc.tsv

    echo ${fields} > datasets/${taxon}_genomes.csv
    cat ${taxon}_genomes*.tsv | awk '{if (\$NF < 15000000 ) print \$0}' | sort | uniq | tr "\\t" "," >> datasets/${taxon}_genomes.csv

    echo "done"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets --version | awk '{print \$NF}' )
        dataformat: \$(dataformat version | awk '{print \$NF}' )
    END_VERSIONS
  """
}

// It is faster if datasets can download the entire list at a time, but there is a 20 minute timeout for downloading.
// The '||' is to allow each genome to be downloaded on its own, which is longer overall but each genome should be less than 20 minutes.
process DATASETS_DOWNLOAD {
  tag           "Downloading Genomes"
  label         "process_medium"
  container     'staphb/ncbi-datasets:16.35.0'
  
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
    done  

    # removing MAGS
    rm -rf genomes/*:_*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      datasets: \$(datasets --version | awk '{print \$NF}')
    END_VERSIONS
  """
}

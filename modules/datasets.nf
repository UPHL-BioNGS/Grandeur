process datasets_summary {
  tag "${taxon}"

  input:
  val(taxon)

  output:
  file("datasets/${taxon}*genomes.txt")                         , emit: genomes
  path "logs/${task.process}/${taxon}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p datasets logs/!{task.process}
    log_file=logs/!{task.process}/!{taxon}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    datasets --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    datasets summary genome taxon "!{taxon}" --reference --as-json-lines | \
        dataformat tsv genome --fields accession,assminfo-refseq-category,organism-name --elide-header | \
        grep representative | \
        head -n !{params.fastani_max_genomes} > \
        datasets/!{taxon}_representative_genomes.txt

    if [ ! -s genome_ids.txt ] 
    then

        echo "representative genome for !{taxon} taxon not found"
        datasets summary genome taxon "!{taxon}" --reference --as-json-lines | \
            dataformat tsv genome --fields accession,assminfo-refseq-category,organism-name --elide-header | \
            head -n !{params.fastani_max_genomes} > \
            datasets/!{taxon}_genomes.txt
    fi

    exit
    '''
}

process datasets_download {
  tag "${taxon}"
  label "medcpus"

  input:
  file(ids)

  output:
  file("datasets/fastani_refs.tar.gz")                                   , emit: tar
  path "logs/${task.process}/datasets_download.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p datasets genomes logs/!{task.process}
    log_file=logs/!{task.process}/datasets_download.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    datasets --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    sort !{ids} | uniq | cut -f 1 > id_list.txt

    datasets download genome accession --inputfile id_list.txt --filename ncbi_datasets.zip

    unzip ncbi_datasets.zip
    
    while read line
    do
        accesion="$(echo $line | cut -f 1 -d , )"
        organism="$(echo $line | cut -f 3 -d , | sed 's/ /_/g')"
        cat ncbi_datasets/data/$accesion/$accesion*.fna | sed 's/ /_/g' | sed 's/,//g' > genomes/${organism}_${accesion}.fna
    done < <(sort !{ids} | uniq | tr "\t" "," )

    tar -czvf datasets/fastani_refs.tar.gz genomes/

    exit
    '''
}
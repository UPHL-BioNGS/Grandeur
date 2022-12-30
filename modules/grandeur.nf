process species {
  tag "Creating list of species"

  input:
  file(results)

  output:
  path "datasets/species_list.txt"                                     , emit: species                                    
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p datasets logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    if [ -f "blobtools_species.txt" ]
    then 
        cut -f 2 blobtools_species.txt | grep -v no-hit | grep -v undef | grep -v name | grep -v ".sp" >> species.txt
    fi

    if [ -f "kraken2_summary.csv" ]
    then 
        cut -f 8 -d , kraken2_summary.csv | grep -v name | grep -v sp. >> species.txt
    fi

    sort species.txt | uniq > datasets/species_list.txt
  '''
}

process decompression {
  tag "Decompressing genome file"

  input:
  file(compressed)

  output:
  path "genomes"                                                       , emit: decompressed                                    
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    tar -xvf !{compressed}
    rm !{compressed}

    if [ ! -d "genomes" ]
    then
      filename=$(ls -d * | grep -v logs)
      mv $filename genomes
    fi
  '''
}
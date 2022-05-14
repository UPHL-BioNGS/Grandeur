process create {
  tag "${sample}"

  input:
  tuple val(sample), file(contig), file(blastn), file(bam)

  output:
  tuple val(sample), file("blobtools/${sample}.blobDB.json")            , emit: file
  path "blobtools/${sample}.${sample}.sorted.bam.cov"                   , emit: files
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools create !{params.blobtools_create_options} \
      -o blobtools/!{sample} \
      -i !{contig} \
      -b !{sample}.sorted.bam \
      -t !{blastn} \
      2>> $err_file >> $log_file
  '''
}

process view {
  tag "${sample}"

  input:
  tuple val(sample), file(json)

  output:
  tuple val(sample), file("blobtools/${sample}.blobDB.table.txt")       , emit: file
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools view !{params.blobtools_view_options} \
      -i !{json} \
      -o blobtools/ \
      2>> $err_file >> $log_file
  '''
}

process blobtools {
  tag "${sample}"

  input:
  tuple val(sample), file(json)

  output:
  path "blobtools/${sample}.*"                                          , emit: files
  tuple val(sample), env(blobtools_species)                             , emit: species
  tuple val(sample), env(blobtools_perc)                                , emit: perc
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p blobtools logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "blobtools version $(blobtools -v)" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    blobtools plot !{params.blobtools_plot_options} \
      -i !{json} \
      -o blobtools/ \
      2>> $err_file 2>> $log_file

    perc='0.0'
    blobtools_species='missing'
    while read line
    do
      new_perc=$(echo $line | cut -f 13 -d " " | sed 's/%//g')
      min=$(echo $perc $new_perc | awk '{if ($1 > $2) print $1; else print $2}')
      if [ "$min" != "$perc" ]
      then
        perc=$new_perc
        blobtools_species=$(echo $line | cut -f 1 -d " " )
        blobtools_perc=$(echo $line | cut -f 13 -d " " )
      fi
    done < <(grep -vw all blobtools/!{sample}*.stats.txt | grep -v "# name" | tr ' ' '_' | grep '%')
  '''
}

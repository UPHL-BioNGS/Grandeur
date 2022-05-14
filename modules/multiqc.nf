process multiqc {
  tag "multiqc"

  when:
  params.fastq_processes =~ /multiqc/ || params.contig_processes =~ /multiqc/

  input:
  file(seqyclean)
  file(kraken2_contigs)
  file(quast)
  file(fastqc)
  file(kraken2_fastq)
  file(prokka)

  output:
  path "${task.process}/multiqc_report.html", optional: true                 , emit: report
  path "${task.process}/multiqc_data/*"     , optional: true                 , emit: data_folder
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}", emit: log_files

  shell:
  '''
    mkdir -p !{task.process} quast fastqc kraken2 prokka seqyclean logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    multiqc --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    for quast_file in !{quast}
    do
      if [ -f "$quast_file" ]
      then
        sample=$(echo $quast_file | sed 's/_quast_report.tsv//g' | head -n 1 )
        mkdir -p quast/$sample
        mv $quast_file quast/$sample/report.tsv
      fi
    done

    multiqc !{params.multiqc_options} \
      --outdir !{task.process} \
      --cl_config "prokka_fn_snames: True"  \
      . \
      2>> $err_file >> $log_file
  '''
}

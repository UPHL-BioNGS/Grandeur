process plasmidfinder {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/plasmidfinder:2.1.6'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3

  input:
  tuple val(sample), file(file)

  output:
  path "plasmidfinder/${sample}/*"                                    , emit: files
  path "plasmidfinder/${sample}_plasmidfinder.tsv"                    , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"     , emit: log

  shell:
  '''
    mkdir -p plasmidfinder/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "plasmidfinder.py: no version" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    plasmidfinder.py !{params.plasmidfinder_options} \
      -i !{file} \
      -o plasmidfinder/!{sample} \
      --extented_output \
      | tee -a $log_file

    head -n 1 plasmidfinder/!{sample}/results_tab.tsv | awk '{print "sample\\t" $0 }' > plasmidfinder/!{sample}_plasmidfinder.tsv
    tail -n +2 plasmidfinder/!{sample}/results_tab.tsv | awk -v sample=!{sample} '{print sample "\\t" $0 }' >> plasmidfinder/!{sample}_plasmidfinder.tsv
  '''
}

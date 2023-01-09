process serotypefinder {
  tag           "${sample}"
  label         "medcpus"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/serotypefinder:2.0.1'
  maxForks      10
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(file), val(flag)

  output:
  path "serotypefinder/${sample}/*"                                     , emit: files
  path "serotypefinder/${sample}_serotypefinder.tsv"                    , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"       , emit: log

  shell:
  '''
    mkdir -p serotypefinder/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    serotypefinder.py !{params.serotypefinder_options} \
      -i !{file} \
      -o serotypefinder/!{sample} \
      -x \
      | tee -a $log_file

    cp serotypefinder/!{sample}/results_tab.tsv serotypefinder/!{sample}_serotypefinder.tsv

    head -n 1 serotypefinder/!{sample}/results_tab.tsv | awk '{print "sample\\t" $0 }' > serotypefinder/!{sample}_serotypefinder.tsv
    tail -n +2 serotypefinder/!{sample}/results_tab.tsv | awk -v sample=!{sample} '{print sample "\\t" $0 }' >> serotypefinder/!{sample}_serotypefinder.tsv
  '''
}

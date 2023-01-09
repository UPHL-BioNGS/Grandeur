process kleborate {
  tag           "${sample}"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/kleborate:2.1.0'
  maxForks      10
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(contig), val(flag)

  output:
  path "kleborate/${sample}_results.csv"                         , emit: collect
  path "kleborate/${sample}_results.txt"                         , emit: result
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p kleborate logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    kleborate --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    kleborate !{params.kleborate_options} \
      -o kleborate/!{sample}_results.txt \
      -a !{contig} \
      | tee -a $log_file

    head -n 1 kleborate/!{sample}_results.txt | tr "\\t" "," | awk '{print "sample," $0}' > kleborate/!{sample}_results.csv
    tail -n 1 kleborate/!{sample}_results.txt | tr "\\t" "," | awk -v sample=!{sample} '{print sample "," $0}' >> kleborate/!{sample}_results.csv
  '''
}

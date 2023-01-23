process kleborate {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/kleborate:2.1.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(contig), val(flag)

  output:
  path "kleborate/${sample}_results.tsv"                         , emit: collect
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

    head -n 1 kleborate/!{sample}_results.txt | awk '{print "sample\\t" $0}' > kleborate/!{sample}_results.tsv
    tail -n 1 kleborate/!{sample}_results.txt | awk -v sample=!{sample} '{print sample "\\t" $0}' >> kleborate/!{sample}_results.tsv
  '''
}

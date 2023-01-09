process shigatyper {
  tag           "${sample}"
  label         "medcpus"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/shigatyper:2.0.1'
  maxForks      10
  
  when:
  flag =~ 'found'

  input:
  tuple val(sample), file(input), val(flag)

  output:
  path "shigatyper/${sample}_shigatyper.tsv"                     , emit: files
  path "shigatyper/${sample}_shigatyper-hits.tsv", optional: true, emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p shigatyper logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    shigatyper --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    shigatyper !{params.shigatyper_options} \
      --SE !{input} \
      --name !{sample} \
      | tee -a $log_file

    hits=$(find shigatyper -iname "*hits.tsv" | head -n 1)
    if [ -f "$hits" ]
    then
      head -n 1  $hits | awk '{print "sample\\t" $0}' > shigatyper/!{sample}_shigatyper-hits.tsv
      tail -n +2 $hits | awk -v sample=!{sample} '{print sample "\\t" $0}' >> shigatyper/!{sample}_shigatyper-hits.tsv
      rm $hits
    fi

    cat *tsv > shigatyper/!{sample}_shigatyper.tsv
  '''
}

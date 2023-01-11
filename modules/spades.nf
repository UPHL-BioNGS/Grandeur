process spades {
  tag           "${sample}"
  label         "maxcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/spades:3.15.5'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  //#UPHLICA cpus   8
  
  input:
  tuple val(sample), file(reads)

  output:
  path "spades/${sample}/*"                                              , emit: files
  tuple val(sample), file("contigs/${sample}_contigs.fa"), optional: true, emit: contigs
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"        , emit: log

  shell:
  '''
    mkdir -p spades contigs logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    spades.py --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    spades.py !{params.spades_options} \
      -1 !{reads[0]} \
      -2 !{reads[1]} \
      --threads !{task.cpus} \
      -o spades/!{sample} \
      | tee -a $log_file

    if [ -f "spades/!{sample}/contigs.fasta" ] ; then cp spades/!{sample}/contigs.fasta contigs/!{sample}_contigs.fa ; fi
  '''
}

process fastani {
  tag           "${sample}"
  label         "medcpus"
  stageInMode   "copy"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/fastani:1.33'
  maxForks      10
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  
  input:
  tuple val(sample), file(contigs), path(genomes)

  output:
  tuple val(sample), file("fastani/${sample}_fastani.csv")       , emit: results
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p fastani logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "fastANI version: " >> $log_file
    fastANI --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file
    
    find !{genomes} -iname "*.fna" -o -iname "*.fasta" -o -iname "*.fa" > reference_list.txt

    fastANI !{params.fastani_options} \
      --threads !{task.cpus} \
      -q !{contigs} \
      --rl reference_list.txt \
      -o fastani/!{sample}.txt \
      | tee -a $log_file

    echo "sample,query,reference,ANI estimate,total query sequence fragments,fragments aligned as orthologous matches" > fastani/!{sample}_fastani.csv
    cat fastani/!{sample}.txt | tr "\\t" "," | awk -v sample=!{sample} '{ print sample "," $0 }' >> fastani/!{sample}_fastani.csv
  '''
}

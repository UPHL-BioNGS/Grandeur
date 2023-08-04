process snp_dists {
  tag "SNP matrix"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/snp-dists:0.8.2'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '24h'
  
  input:
  tuple file(contigs), val(num_samples), val(num_core_genes)

  output:
  path "snp-dists/snp_matrix.txt"                                      , emit: snp_matrix
  path "snp-dists/snp_matrix_with_qc.txt"
  path "snp-dists/roary_metrics_mqc.csv"                               , emit: for_multiqc                                      
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p snp-dists logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    snp-dists -v >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    snp-dists !{params.snp_dists_options} \
      !{contigs} \
      > snp-dists/snp_matrix.txt

    genome_length=$(cat !{contigs} | tr "\n" ";" | sed 's/>[^>]*//2g' | tr ";" "\n" | grep -v ">" | wc -c )

    sed '0,/,/s/,/num_samples=!{num_samples};num_core_genes=!{num_core_genes},/' snp-dists/snp_matrix.txt > snp-dists/snp_matrix_with_qc.txt

    echo "num_samples,num_core_genes,core_genome_length"     >  snp-dists/roary_metrics_mqc.csv
    echo "!{num_samples},!{num_core_genes},${genome_length}" >> snp-dists/roary_metrics_mqc.csv
  '''
}

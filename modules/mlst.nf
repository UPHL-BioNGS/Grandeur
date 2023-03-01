process mlst {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mlst:2.23.0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  tuple val(sample), file(contig)

  output:
  path "mlst/${sample}_mlst.tsv"                                 , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p mlst logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    mlst --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    echo -e "sample\\tfilename\\tmatching PubMLST scheme\\tST\\tID1\\tID2\\tID3\\tID4\\tID5\\tID6\\tID7\\tID8\\tID9\\tID10\\tID11\\tID12\\tID13\\tID14\\tID15" > mlst/!{sample}_mlst.tsv

    mlst !{params.mlst_options} \
      --threads !{task.cpus} \
      !{contig} | \
      awk -v sample=!{sample} '{print sample "\\t" $1 "\\t" $2 "\\t" $3 "\\t" $4 "\\t" $5 "\\t" $6 "\\t" $7 "\\t" $8 "\\t" $9 "\\t" $10 "\\t" $11 "\\t" $12 "\\t" $13 "\\t" $14 "\\t" $15 "\\t" $16 "\\t" $17 "\\t" $18}' >> mlst/!{sample}_mlst.tsv
  '''
}

process mlst {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mlst:2.22.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3

  input:
  tuple val(sample), file(contig)

  output:
  path "mlst/${sample}_mlst.csv"                                 , emit: collect
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

    echo "sample,filename,matching PubMLST scheme,ST,ID1,ID2,ID3,ID4,ID5,ID6,ID7,ID8,ID9,ID10,ID11,ID12,ID13,ID14,ID15" > mlst/!{sample}_mlst.csv

    mlst !{params.mlst_options} \
      --threads !{task.cpus} \
      !{contig} | \
      awk -v sample=!{sample} '{print sample "," $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7 "," $8 "," $9 "," $10 "," $11 "," $12 "," $13 "," $14 "," $15 "," $16 "," $17 "," $18}' >> mlst/!{sample}_mlst.csv
  '''
}

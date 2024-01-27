process fastani {
  tag           "$meta.id"
  label         "process_medium"
  stageInMode   "copy"
  publishDir    path: params.outdir, mode: 'copy', pattern: 'logs/*/*log'
  publishDir    path: params.outdir, mode: 'copy', pattern: 'fastani/*' 
  container     'staphb/fastani:1.34'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  //#UPHLICA memory 26.GB
  //#UPHLICA cpus 7
  //#UPHLICA time '10m'
  
  input:
  tuple val(meta), file(contigs), file(genomes)

  output:
  tuple val(meta), file("fastani/${prefix}_fastani.csv")       , emit: results, optional: true
  tuple val(meta), env(top_hit), path("top_hit/*")             , emit: top_hit, optional: true
  path "fastani/*"                                               , emit: everything
  path "logs/${task.process}/${prefix}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p fastani logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : ${task.container}" >> $log_file
    echo "fastANI version: " >> $log_file
    fastANI --version 2>> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    echo ${genomes} | tr " " "\\n" | sort > reference_list.txt

    fastANI ${params.fastani_options} \
      --threads ${task.cpus} \
      -q ${contigs} \
      --rl reference_list.txt \
      -o fastani/${prefix}.txt \
      | tee -a $log_file

    echo "sample,query,reference,ANI estimate,total query sequence fragments,fragments aligned as orthologous matches" > fastani/${prefix}_fastani.csv
    cat fastani/${prefix}.txt | sed 's/,//g' | tr "\\t" "," | awk -v sample=${prefix} '{ print sample "," $0 }' >> fastani/${prefix}_fastani.csv

    top_hit=$(head -n 2 fastani/${prefix}_fastani.csv | tail -n 1 | cut -f 3 -d , )
    if [ -f "$top_hit" ]
    then
      mkdir -p top_hit
      cp $top_hit top_hit/.
      gzip -d top_hit/*.gz || ls top_hit
      chmod 664 top_hit/*
    fi
  """
}

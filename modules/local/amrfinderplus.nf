process amrfinderplus {
  tag           "$meta.id"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/ncbi-amrfinderplus:3.11.26-2023-11-15.1'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA cpus 7
  //#UPHLICA memory 26.GB
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  //#UPHLICA time '10m'

  input:
  tuple val(sample), file(contigs), val(genus), val(species)

  output:
  path "ncbi-AMRFinderplus/${sample}_amrfinder_plus.txt"         , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log
  path  "versions.yml"                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
      def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
  '''
    mkdir -p ncbi-AMRFinderplus logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file 
    amrfinder --version >> $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    organism=$(amrfinder -l | tr " " "\\n" | grep -i !{genus} | grep -i !{species} | sed 's/,//g' | head -n 1 )
    if [ -z "$organism" ] ; then organism=$(amrfinder -l | tr " " "\\n" | grep -i !{genus} | sed 's/,//g' | head -n 1 ) ; fi
    if [ -n "$organism" ]
    then
      organism_check="--organism $organism"
      echo "Top organism result of !{genus} !{species} matched with $organism" >> $log_file
    elif [ "!{genus}" == "Shigella" ]
    then
      organism_check="--organism Escherichia"
      echo "--organism Escherichia with be used because of top organism result of !{genus}" >> $log_file
    else
      organism_check=''
      echo "Top organism result of !{genus} !{species} did not match any of the organisms" >> $log_file
    fi

    amrfinder !{params.amrfinderplus_options} \
      --nucleotide !{contigs} \
      --threads !{task.cpus} \
      --name !{sample} \
      --output ncbi-AMRFinderplus/!{sample}_amrfinder_plus.txt \
      $organism_check \
      --plus \
      | tee -a $log_file
  '''
}

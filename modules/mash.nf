process mash {
  tag "${sample}"
  label "medcpus"

  input:
  tuple val(sample), file(reads)

  output:
  tuple val(sample), env(genome_size)                                  , emit: genome_size
  tuple val(sample), env(coverage)                                     , emit: coverage
  tuple val(sample), file("mash/${sample}_mashdist.txt"), optional: true, emit: files
  tuple val(sample), env(genus)                                         , emit: genus
  tuple val(sample), env(species)                                       , emit: species
  tuple val(sample), env(full_mash)                                     , emit: full
  tuple val(sample), env(pvalue)                                        , emit: pvalue
  tuple val(sample), env(distance)                                      , emit: distance
  tuple val(sample), env(salmonella_flag)                               , emit: salmonella_flag
  tuple val(sample), env(ecoli_flag)                                    , emit: ecoli_flag
  tuple val(sample), env(klebsiella_flag)                               , emit: klebsiella_flag
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log,err"       , emit: log


  shell:
  '''
    mkdir -p mash logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "mash version: $(mash --version | head -n 1 )" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    cat !{reads} | mash sketch !{params.mash_sketch_options} -o mash/!{sample} - 2>> $err_file | tee -a $log_file



    genome_size=$(grep "Estimated genome size" $err_file | awk '{print $4}' )
    coverage=$(grep "Estimated coverage" $err_file | awk '{print $3}' )

    if [ -z "$genome_size" ] ; then genome_size=0 ; fi
    if [ -z "$coverage" ] ; then coverage=0 ; fi


    mash dist -p !{task.cpus} !{params.mash_options} !{params.mash_reference} mash/!{sample}.msh | sort -gk3 > mash/!{sample}_mashdist.txt

    if [ ! -s "mash/!{sample}_mashdist.txt" ]
    then
      echo "!{sample} had no mash results with '!{params.mash_options}'. Trying again without those parameters." | tee -a $log_file
      mash dist -p !{task.cpus} !{params.mash_reference} !{msh} | sort -gk3 > mash/!{sample}_mashdist.txt
    fi

    mash_result=($(head -n 1 mash/!{sample}_mashdist.txt | head -n 1 | cut -f 1 | cut -f 8 -d "-" | cut -f 1,2 -d "_" | cut -f 1 -d "." | tr "_" " " ) 'missing' 'missing')
    genus=${mash_result[0]}
    species=${mash_result[1]}
    full_mash=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 1 )
    pvalue=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 4 )
    distance=$(head -n 1 mash/!{sample}_mashdist.txt | cut -f 3 )
    if [ -z "$full_mash" ] ; then full_mash='missing' ; fi
    if [ -z "$pvalue" ] ; then pvalue='NA' ; fi
    if [ -z "$distance" ] ; then distance='NA' ; fi

    salmonella_flag=''
    find_salmonella=$(head mash/!{sample}_mashdist.txt | grep "Salmonella" | head -n 1 )
    if [ -n "$find_salmonella" ] ; then salmonella_flag="found" ; else salmonella_flag="not" ; fi

    ecoli_flag=''
    find_ecoli=$(head mash/!{sample}_mashdist.txt | grep -e "Escherichia" -e "Shigella" | head -n 1 )
    if [ -n "$find_ecoli" ] ; then ecoli_flag="found" ; else ecoli_flag="not" ; fi

    klebsiella_flag=''
    find_klebsiella=$(head mash/!{sample}_mashdist.txt | grep -e "Klebsiella" -e "Enterobacter" -e "Serratia" | head -n 1 )
    if [ -n "$find_klebsiella" ] ; then klebsiella_flag="found" ; else klebsiella_flag="not" ; fi
  '''
}

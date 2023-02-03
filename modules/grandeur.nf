process species {
  tag           "Creating list of species"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.1.5'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'
  
  input:
  file(results)

  output:
  path "datasets/species_list.txt"                                     , emit: species                                    
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p datasets logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    if [ -f "blobtools_species.txt" ]
    then 
      cut -f 2 blobtools_species.txt >> species.txt
    fi

    if [ -f "kraken2_summary.csv" ]
    then 
      cut -f 8 -d , kraken2_summary.csv >> species.txt
    fi

    if [ -f "mash_summary.csv" ]
    then
      cut -f 7 -d , mash_summary.csv >> species.txt
    fi

    grep -v no-hit species.txt | grep -v undef | grep -v name | grep -v "sp." | sort | uniq > datasets/species_list.txt
  '''
}

process decompression {
  tag           "Decompressing genome file"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.1.5'
  maxForks      10
  stageInMode   'copy'
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  file(compressed)

  output:
  path "genomes"                                                       , emit: decompressed                                    
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    zcat !{compressed} | tar -xvf -
    rm !{compressed}

    if [ ! -d "genomes" ]
    then
      filename=$(ls -d * | grep -v logs)
      mv $filename genomes
    fi
  '''
}

process flag {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.1.5'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  tuple val(sample), file(files)

  output:
  tuple val(sample), env(salmonella_flag)                        , emit: salmonella_flag
  tuple val(sample), env(klebsiella_flag)                        , emit: klebsiella_flag
  tuple val(sample), env(ecoli_flag)                             , emit: ecoli_flag
  tuple val(sample), env(genus), env(species)                    , emit: organism
  path "flag/${sample}_flag.csv"                                 , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p flag logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    if [ -f "!{sample}_fastani.csv" ]
    then 
      awk -F "," '{if ($4 > 90) print $0}' !{sample}_fastani.csv > smaller_fastani.csv
      genus=$(head   -n 2 !{sample}_fastani.csv | tail -n 1 | cut -f 2 -d , | cut -f 2 -d '/' | cut -f 1 -d "_")
      species=$(head -n 2 !{sample}_fastani.csv | tail -n 1 | cut -f 2 -d , | cut -f 2 -d '/' | cut -f 2 -d "_")
    else
      touch smaller_fastani.csv
      genus="unknown"
      species="unknown"
    fi

    files=$(ls !{files} | grep -v fastani)

    echo "Looking for Salmonella:" >> $log_file
    salmonella_flag=''
    find_salmonella=$(head -n 10 $files smaller_fastani.csv | grep "Salmonella" | tee -a $log_file | head -n 1 )
    if [ -n "$find_salmonella" ] ; then salmonella_flag="found" ; else salmonella_flag="not" ; fi
    
    echo "Looking for E. coli and Shigella:" >> $log_file
    ecoli_flag=''
    find_ecoli=$(head -n 10 $files smaller_fastani.csv | grep -e "Escherichia" -e "Shigella" | tee -a $log_file | head -n 1 )
    if [ -n "$find_ecoli" ] ; then ecoli_flag="found" ; else ecoli_flag="not" ; fi

    echo "Looking for Klebsiella:" >> $log_file
    klebsiella_flag=''
    find_klebsiella=$(head -n 10 $files smaller_fastani.csv | grep -e "Klebsiella" -e "Enterobacter" -e "Serratia" | tee -a $log_file | head -n 1 )
    if [ -n "$find_klebsiella" ] ; then klebsiella_flag="found" ; else klebsiella_flag="not" ; fi

    if [ -z "$genus" ]   ; then genus=unknown ; fi
    if [ -z "$species" ] ; then species=unknown ; fi

    echo "sample,genus,species,salmonella_flag,ecoli_flag,klebsiella_flag" > flag/!{sample}_flag.csv
    echo "!{sample},$genus,$species,$salmonella_flag,$ecoli_flag,$klebsiella_flag" >> flag/!{sample}_flag.csv
  '''
}

process size {
  tag           "${sample}"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.1.5'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'
    
  input:
  tuple val(sample), file(mash_err), file(fastani), val(top_hit), file(genome_sizes), file(datasets_summary)

  output:
  tuple val(sample), env(size)                                   , emit: size
  tuple val(sample), env(genus), env(species), env(accession)    , emit: organism
  path "size/${sample}_size.csv"                                 , emit: collect
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p size logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    echo "This process is still here for a placeholder. Sometimes when contamination occurs, the genome size from mash is very different than the expectation. TBC" >> $log_file

    mash_size="NA"
    expected_size="NA"
    top_hit_size="NA"
    ncbi_size="NA"
    size=''

    genus=$(echo !{top_hit} | cut -f 2 -d "/" | cut -f 1 -d "_" )
    species=$(echo !{top_hit} | cut -f 2 -d "_" )
    accession=$(echo !{top_hit} | sed 's/.*_GC/GC/g' | cut -f 1,2 -d '.' )

    top_hit_size=$(grep -v ">" !{fastani} | wc -c )
    size=$top_hit_size

    if [ -f "datasets_summary.csv" ]
    then
      ncbi_size=$(grep $accession datasets_summary.csv | cut -f 5 -d "," | head -n 1 )
      echo "The expected size based on the fastANI top hit is $ncbi_size" | tee -a $log_file
      size=$ncbi_size
    fi

    genus_check=$(grep -v '#' !{genome_sizes} | grep $genus | head -n 1)
    if [ -n "$genus_check" ]
    then
      species_check=$(grep -v '#' !{genome_sizes} | grep $genus | grep $species | head -n 1 )
      
      if [ -n "$species_check" ]
      then
        expected_size=$(grep -v '#' !{genome_sizes} | grep $genus | grep $species | awk '{print $3}' | sed 's/,$//g' | head -n 1 )
        echo "The expected size based on $genus and $species is $expected_size" | tee -a $log_file
        size=$expected_size
      else
        expected_size=$(grep -v '#' !{genome_sizes} | grep $genus | awk '{print $3}' | sed 's/,$//g' | head -n 1 )
        echo "The expected size based on using $genus is $expected_size" | tee -a $log_file
        size=$expected_size
      fi
    else
      expected_size="not found"
      echo "The expected size based on $genus and $species was not found" | tee -a $log_file
    fi

    mash_size="$(grep "Estimated genome size" !{mash_err} | awk '{print $4 }' | sort -gr | tr '\\n' ' ' )"
    echo "The expected size based on kmers from mash is $mash_size" | tee -a $log_file
    
    echo "sample,genus,species,accession,size,expected,top_hit,ncbi,mash" > size/!{sample}_size.csv
    echo "!{sample},$genus,$species,$accession,$size,$expected_size,$top_hit_size,$ncbi_size,$mash_size" >> size/!{sample}_size.csv

    if [ -z "$size" ] ; then size=$(echo $mash_size | cut -f 1 -d " " ) ; fi
    
    echo "The final size is $size" | tee -a $log_file
  '''
}

process representative {
  tag           "${accession}"
  container     'quay.io/biocontainers/pandas:1.1.5'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  tuple val(accession), path(genomes)

  output:
  tuple path("representative/*"), env(genus), env(species), val(accession), emit: representative
  path "logs/${task.process}/${accession}.${workflow.sessionId}.log"      , emit: log_files

  shell:
  '''
    mkdir -p representative logs/!{task.process}
    log_file=logs/!{task.process}/!{accession}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    fasta=$(ls !{genomes}/*!{accession}* | head -n 1 | cut -f 2 -d "/")
    genus=$(echo $fasta | cut -f 1 -d "_" )
    species=$(echo $fasta | cut -f 2 -d "_" )

    cp !{genomes}/*!{accession}* representative/.
  '''
}
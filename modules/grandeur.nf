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
      cut -f 7 -d , mash_summary.csv | tail -n+2 >> species.txt
    fi

    grep -v no-hit species.txt | grep -v undef | grep -v name | sort | uniq > datasets/species_list.txt
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
  tuple val(sample), env(streppneu_flag)                         , emit: streppneu_flag
  tuple val(sample), env(legionella_flag)                        , emit: legionella_flag
  tuple val(sample), env(klebacin_flag)                          , emit: klebacin_flag
  tuple val(sample), env(strepa_flag)                            , emit: strepa_flag
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

    
    echo "Looking for Strep A organisms:" >> $log_file
    strepa_flag=''
    find_strepa=$(head -n 10 $files smaller_fastani.csv | grep "Streptococcus" | grep -e "pyogenes" -e "dysgalactiae" -e "anginosus" | tee -a $log_file | head -n 1 )
    if [ -n "$find_strepa" ] ; then strepa_flag='found' ; else strepa_flag='not' ; fi

    echo "Looking for Streptococcus pneumoniae organisms:" >> $log_file
    streppneu_flag_flag=''
    find_streppneu=$(head -n 10 $files smaller_fastani.csv | grep "Streptococcus" | grep "pneumoniae" | tee -a $log_file | head -n 1 )
    if [ -n "$find_streppneu" ] ; then streppneu_flag='found' ; else streppneu_flag='not' ; fi
  
    echo "Looking for Legionella organisms:" >> $log_file
    legionella_flag=''
    find_legionella=$(head -n 10 $files smaller_fastani.csv | grep "Legionella" | tee -a $log_file | head -n 1 )
    if [ -n "$find_legionella" ] ; then legionella_flag='found' ; else legionella_flag='not' ; fi

    echo "Looking for Klebsiella or Acinetobacter:" >> $log_file
    klebacin_flag=''
    if [ -n "$find_klebsiella" ]
    then
      klebacin_flag='found'
    else
      find_acin=$(head -n 10 $files smaller_fastani.csv | grep "Acinetobacter" | tee -a $log_file | head -n 1 )
      if [ -n "$find_acin" ] ; then klebacin_flag='found' ; else klebacin_flag='not' ; fi
    fi

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
  tuple val(sample), file(results)
  // results should include the mash results, the mash err files, the fastani results, the genome sizes file, the genomes file, the quast file

  output:
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

    # Step 1 : Getting the top hit
    genus=""
    species=""
    accession=""

    if [ -f "!{sample}_fastani.csv" ] 
    then
      if [ "$(wc -l !{sample}_fastani.csv | awk '{print $1}' )" -gt 1 ]
      then
        genus=$(head     -n 2 !{sample}_fastani.csv | tail -n 1 | cut -f 3 -d "," | cut -f 1 -d "_" )
        species=$(head   -n 2 !{sample}_fastani.csv | tail -n 1 | cut -f 3 -d "," | cut -f 2 -d "_" )
        accession=$(head -n 2 !{sample}_fastani.csv | tail -n 1 | sed 's/.*_GC/GC/g' | cut -f 1,2 -d '.' )
      fi
    fi

    if [ -z "$genus" ] && [ -f "!{sample}.summary.mash.csv" ]
    then
      genus=$(head -n 2 !{sample}.summary.mash.csv | tail -n 1 | cut -f 7 -d "," | cut -f 1 -d "_" )
      species=$(head -n 2 !{sample}.summary.mash.csv | tail -n 1 | cut -f 7 -d "," | cut -f 2 -d "_" )
      accession=$(head -n 2 !{sample}.summary.mash.csv | tail -n 1 | cut -f 2 -d "," | cut -f 1,2 -d '_' )
    fi
    
    # Step 2 : Using this information to get the estimated genome size
    datasets_size=""
    expected_size=""
    mash_size=""
    quast_size=""

    if [ -f "datasets_summary.csv" ] && [ -n "$accession" ]
    then
      datasets_size=$(grep $accession datasets_summary.csv | cut -f 5 -d "," | head -n 1 )
      if [ -z "$datasets_size" ] && [ -n "$species" ] ; then datasets_size=$(grep $genus datasets_summary.csv | grep $species | cut -f 5 -d "," | head -n 1 ) ; fi
      echo "The expected size based on the fastANI top hit is $datasets_size" | tee -a $log_file
    fi

    if [ -n "$genus" ]
    then
      genus_check=$(grep -v '#' genome_sizes.json | grep $genus | head -n 1)
      if [ -n "$genus_check" ]
      then
        species_check=$(grep -v '#' genome_sizes.json | grep $genus | grep $species | head -n 1 )
      
        if [ -n "$species_check" ]
        then
          expected_size=$(grep -v '#' genome_sizes.json | grep $genus | grep $species | awk '{print $3}' | sed 's/,$//g' | head -n 1 )
          echo "The expected genome size based on $genus and $species is $expected_size" | tee -a $log_file
        else
          expected_size=$(grep -v '#' genome_sizes.json | grep $genus | awk '{print $3}' | sed 's/,$//g' | head -n 1 )
          echo "The expected genome size based on using $genus is $expected_size" | tee -a $log_file
        fi
      else
        echo "The expected genome size based on $genus and $species was not found" | tee -a $log_file
      fi
    else
      echo "There is no genus to look for" | tee -a $log_file
    fi

    mash_header="mash_size"
    if [ -f "!{sample}.!{workflow.sessionId}.err" ]
    then
      mash_check=$(grep "Estimated genome size" !{sample}.!{workflow.sessionId}.err | head -n 1)
      
      if [ -n "$mash_check" ]
      then
        mash_sizes=($(grep "Estimated genome size" !{sample}.!{workflow.sessionId}.err | awk '{print $4 }'))
        i=1
      
        for err_size in ${mash_sizes[@]}
        do
          err_size=$(printf "%.0f" $err_size)
          if [ -z "$mash_size" ]
          then
            mash_size="$err_size"
          else
            mash_size="$mash_size,$err_size"
            mash_header="$mash_header,mash_size$i"
          fi
          i=$((i + 1))
        done
        echo "The expected size based on kmers from mash is ${mash_sizes[@]}" | tee -a $log_file
      fi
    fi

    if [ -f "!{sample}_quast_report.tsv" ]
    then
      quast_size=$(grep "Total length" !{sample}_quast_report.tsv | grep -v "=" | awk '{print $3}')
      echo "The total length from quast is $quast_size" | tee -a $log_file
    fi

    # Step 3. Settling on a final genome size for coverage
    if [ -n "$expected_size" ]
    then
      echo "Using size from datasets : $expected_size" | tee -a $log_file
      size=$expected_size
    elif [ -z "$expected_size" ] && [ -n "$datasets_size" ]
    then
      echo "Using size from genomes file : $datasets_size" | tee -a $log_file
      size=$datasets_size
    elif [ -n "$quast_size" ] && [ -z "$datasets_size" ] && [ -z "$expected_size" ]
    then
      echo "Using size from quast : $quast_size" | tee -a $log_file
      size=$quast_size
    else
      echo "A size could not be determined. Setting to 5M" | tee -a $log_file
      size=5000000
    fi
    echo "The final size is $size" | tee -a $log_file
    
    # Step 4. Putting it in a file
    echo "sample,genus,species,accession,size,datasets_size,expected_size,$mash_header,quast_size" > size/!{sample}_size.csv
    echo "!{sample},$genus,$species,$accession,$size,$datasets_size,$expected_size,$mash_size,$quast_size" >> size/!{sample}_size.csv
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
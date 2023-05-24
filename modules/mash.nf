process mash {
  tag           "${sample}"
  label         "medcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mash:2.3'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize' , value: 'himem-medium'
  //#UPHLICA cpus 8
  //#UPHLICA memory 100.GB
  //#UPHLICA time '10m'

  input:
  tuple val(sample), file(fasta), file(fastq), file(reference)

  output:
  path "mash/${sample}.mashdist.txt"                                 , emit: mashdist
  tuple val(sample), file("mash/${sample}.summary.mash.csv")         , emit: results
  tuple val(sample), file("mash/${sample}.${workflow.sessionId}.err"), optional: true, emit: err
  path "logs/${task.process}/${sample}.${workflow.sessionId}.log"    , emit: log

  shell:
  '''
  mkdir -p mash logs/!{task.process}
  log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
  err_file=mash/!{sample}.!{workflow.sessionId}.err

  # time stamp + capturing tool versions
  date > $log_file
  echo "container : !{task.container}" >> $log_file
  echo "mash version: $(mash --version | head -n 1 )" >> $log_file
  echo "Nextflow command : " >> $log_file
  cat .command.sh >> $log_file
 
  if [[ "!{reference}" != *"input"* ]]
  then 
    reference=!{reference}
  else
    reference=/db/RefSeqSketchesDefaults.msh
  fi

  for file in !{fasta} !{fastq}
  do
    if [[ "$file" == *"fastp"* ]] 
    then
      cat $file >> sample.fastq.gz
      if [ ! -s "sample.fastq.gz" ] ; then rm sample.fastq.gz ; fi
    elif [[ "$file" != *"input"* ]]
    then
      cat $file > sample.fasta
    fi 
  done

  if [ -f "sample.fastq.gz" ]
  then
    echo "Mash on fastq files" | tee -a $log_file
    cat sample.fastq.gz | mash sketch !{params.mash_sketch_options} -o !{sample}.fastq - 2>> $err_file | tee -a $log_file

    mash dist -p !{task.cpus} !{params.mash_dist_options} $reference !{sample}.fastq.msh > mash/!{sample}.mashdist.txt.tmp
  fi

  if [ -f "sample.fasta" ]
  then
    echo "Mash on fasta files" | tee -a $log_file
    mash sketch -o !{sample}.fasta sample.fasta 2>> $err_file | tee -a $log_file

    mash dist -p !{task.cpus} !{params.mash_dist_options} $reference !{sample}.fasta.msh >> mash/!{sample}.mashdist.txt.tmp
  fi

  touch mash/!{sample}.mashdist.txt.tmp

  if [ ! -s "mash/!{sample}.mashdist.txt.tmp" ]
  then
    echo "!{sample} had no mash results with '!{params.mash_dist_options}'. Trying again without those parameters." | tee -a $log_file
        
    if [ -s "!{sample}.fastq.msh" ] ; then mash dist -p !{task.cpus} $reference !{sample}.fastq.msh >> mash/!{sample}.mashdist.txt.tmp ; fi

    if [ -s "!{sample}.fasta.msh" ] ; then mash dist -p !{task.cpus} $reference !{sample}.fasta.msh >> mash/!{sample}.mashdist.txt.tmp ; fi
  fi

  sort -gk3 mash/!{sample}.mashdist.txt.tmp | head -n !{params.mash_max_hits} > mash/!{sample}.mashdist.txt

  echo "sample,reference,query,mash-distance,P-value,matching-hashes,organism" > mash/!{sample}.summary.mash.csv

  if [[ "!{reference}" != *"input"* ]]
  then
    while read line
    do 
      organism=$(echo $line | cut -f 3,4 -d "_" | cut -f 1 -d ".")
      echo $line | sed 's/,//g' | awk -v sample=!{sample} -v org=$organism '{print sample "," $1 "," $2 "," $3 "," $4 "," $5 "," org}' >> mash/!{sample}.summary.mash.csv
    done < mash/!{sample}.mashdist.txt
  else
    while read line
    do
      organism=$(echo $line | cut -f 8 -d "-" | cut -f 1,2 -d "_" | cut -f 1 -d ".")
      echo $line | sed 's/,//g' | awk -v sample=!{sample} -v org=$organism '{print sample "," $1 "," $2 "," $3 "," $4 "," $5 "," org}' >> mash/!{sample}.summary.mash.csv
    done < mash/!{sample}.mashdist.txt
  fi
  '''
}

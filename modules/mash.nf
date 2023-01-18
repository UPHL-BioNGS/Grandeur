process mash {
  tag           "${sample}"
  label         "medcpus"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/mash:2.3'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-large'
  //#UPHLICA cpus   8
  
  input:
  tuple val(sample), file(input), file(reference)

  output:
  path "mash/${input[0]}.mashdist.txt"                                 , emit: mashdist
  tuple val(sample), file("mash/${input[0]}.summary.mash.csv")         , emit: results
  tuple val(sample), file("mash/${input[0]}.${workflow.sessionId}.err"), emit: err
  path "logs/${task.process}/${input[0]}.${workflow.sessionId}.log"    , emit: log

  shell:
  if (params.mash_db) {
    '''
      mkdir -p mash logs/!{task.process}
      log_file=logs/!{task.process}/!{input[0]}.!{workflow.sessionId}.log
      err_file=mash/!{input[0]}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date > $log_file
      echo "container : !{task.container}" >> $log_file
      echo "mash version: $(mash --version | head -n 1 )" >> $log_file
      echo "Nextflow command : " >> $log_file
      cat .command.sh >> $log_file

      cat !{input} | mash sketch !{params.mash_sketch_options} -o !{input[0]} - 2>> $err_file | tee -a $log_file

      mash dist -p !{task.cpus} !{params.mash_dist_options} !{reference} !{input[0]}.msh | sort -gk3 | head -n !{params.mash_max_hits} > mash/!{input[0]}.mashdist.txt

      if [ ! -s "mash/!{input[0]}.mashdist.txt" ]
      then
        echo "!{sample} had no mash results with '!{params.mash_dist_options}'. Trying again without those parameters." | tee -a $log_file
        mash dist -p !{task.cpus} !{reference} !{input[0]}.msh | sort -gk3 | head -n !{params.mash_max_hits} > mash/!{input[0]}.mashdist.txt
      fi

      echo "sample,reference,query,mash-distance,P-value,matching-hashes,organism" > mash/!{input[0]}.summary.mash.csv

      while read line
      do
        organism=$(echo $line | cut -f 3,4 -d "_")
        echo $line | sed 's/,//g' | awk -v sample=!{sample} -v org=$organism '{print sample "," $1 "," $2 "," $3 "," $4 "," $5 "," org}' >> mash/!{input[0]}.summary.mash.csv
      done < mash/!{input[0]}.mashdist.txt
    '''
  } else {
    '''
      mkdir -p mash logs/!{task.process}
      log_file=logs/!{task.process}/!{input[0]}.!{workflow.sessionId}.log
      err_file=mash/!{input[0]}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date > $log_file
      echo "container : !{task.container}" >> $log_file
      echo "mash version: $(mash --version | head -n 1 )" >> $log_file
      echo "Nextflow command : " >> $log_file
      cat .command.sh >> $log_file

      cat !{input} | mash sketch !{params.mash_sketch_options} -o !{input[0]} - 2>> $err_file | tee -a $log_file

      mash dist -p !{task.cpus} !{params.mash_dist_options} /db/RefSeqSketchesDefaults.msh !{input[0]}.msh | sort -gk3 | head -n !{params.mash_max_hits} > mash/!{input[0]}.mashdist.txt

      if [ ! -s "mash/!{input[0]}.mashdist.txt" ]
      then
        echo "!{sample} had no mash results with '!{params.mash_dist_options}'. Trying again without those parameters." | tee -a $log_file
        mash dist -p !{task.cpus} /db/RefSeqSketchesDefaults.msh !{input[0]}.msh | sort -gk3 | head -n !{params.mash_max_hits} > mash/!{input[0]}.mashdist.txt
      fi

      echo "sample,reference,query,mash-distance,P-value,matching-hashes,organism" > mash/!{input[0]}.summary.mash.csv

      while read line
      do
        organism=$(echo $line | cut -f 8 -d "-" | cut -f 1,2 -d "_" | cut -f 1 -d ".")
        echo $line | sed 's/,//g' | awk -v sample=!{sample} -v org=$organism '{print sample "," $1 "," $2 "," $3 "," $4 "," $5 "," org}' >> mash/!{input[0]}.summary.mash.csv
      done < mash/!{input[0]}.mashdist.txt
    '''
  }
}
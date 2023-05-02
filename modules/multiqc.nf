process multiqc {
  tag           "multiqc"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
  maxForks      10
  //#UPHLICA errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  //#UPHLICA pod annotation: 'scheduler.illumina.com/presetSize', value: 'standard-medium'
  //#UPHLICA memory 1.GB
  //#UPHLICA cpus 3
  //#UPHLICA time '10m'

  input:
  file(input)

  output:
  path "multiqc/multiqc_report.html", optional: true                   , emit: report
  path "multiqc/multiqc_data/*"     , optional: true                   , emit: data_folder
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log", emit: log_files

  shell:
  '''
    mkdir -p multiqc quast logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log

    # time stamp + capturing tool versions
    date > $log_file
    echo "container : !{task.container}" >> $log_file
    multiqc --version >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    for quast_file in $(ls *_quast_report.tsv)
    do
      sample=$(echo $quast_file | sed 's/_quast_report.tsv//g' | head -n 1 )
      mkdir -p quast/$sample
      mv $quast_file quast/$sample/report.tsv
    done

    if [ -f 'mash_summary.csv' ]
    then 
      organisms=($(cut -f 7 -d , mash_summary.csv | sort | uniq | grep -v organism ))
      samples=($(cut -f 1 -d , mash_summary.csv | sort | uniq | grep -v sample ))

      echo ${organisms[@]} | tr ' ' ',' | awk '{print "sample," $0}' | tee mash_mqc.csv mash_fasta_mqc.csv

      for sample in ${samples[@]}
      do
        line="$sample"
        fasta_line="$sample"

        for organism in ${organisms[@]}
        do
          num=$(grep ^"$sample," mash_summary.csv | grep ",$organism"$ | grep -v sample.fasta | wc -l | awk '{print $1}' )
          if [ -z "$num" ] ; then num=0 ; fi
          line="$line,$num"

          num=$(grep ^"$sample," mash_summary.csv | grep ",$organism"$ | grep sample.fasta | wc -l | awk '{print $1}' )
          if [ -z "$num" ] ; then num=0 ; fi
          fasta_line="$fasta_line,$num"
        done
        echo $line >> mash_mqc.csv
        echo $fasta_line >> mash_fasta_mqc.csv
      done
    fi

    if [ -f 'blobtools_summary.txt' ]
    then
      organisms=($(cut -f 2 blobtools_summary.txt | grep -v all  | grep -v name | sort | uniq ))
      samples=($(cut -f 1 blobtools_summary.txt | grep -v all | grep -v sample | sort | uniq ))

      echo ${organisms[@]} | tr ' ' '\t' | awk '{print "sample\t" $0}' > blobtools_mqc.tsv

      for sample in ${samples[@]}
      do
        line="$sample"

        for organism in ${organisms[@]}
        do
          num=$(grep -w ^"$sample" blobtools_summary.txt | grep -w "$organism" | cut -f 13 )
          if [ -z "$num" ] ; then num=0 ; fi
          line="$line\t$num"
        done
        echo -e $line | sed 's/,//g' >>  blobtools_mqc.tsv
      done
    fi
 
    if [ -f 'fastani_summary.csv' ] ; then awk -F ',' '{print $1 "_" NR "," $3 "," $4 "," $5 "," $6 }' fastani_summary.csv > fastani_mqc.csv ; fi
    if [ -f 'shigatyper_results.txt' ] ; then awk '{print $1 "_" $2 "\t" $3}' shigatyper_results.txt > shigatyper_mqc.tsv ; fi
    if [ -f 'amrfinderplus.txt' ] ; then sed 's/ /_/g' amrfinderplus.txt | awk '{print $1 "_" NR "\t" $7 "\t" $11 "\t" $12 "\t" $13 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $18 "\t" $20 }'   >  amrfinderplus_mqc.txt  ; fi
    if [ -f 'fastqscan_summary.csv' ] ; then awk -F , '{if ($4 > 40) print $0 ",PASS" ; else if ($4 > 20) print $0 ",TBD" ; else print $0 ",FAIL" }' fastqscan_summary.csv > fastqscan_mqc.csv ; fi
    if [ -f 'kleborate_results.tsv' ] ; then cut -f 1,3-12 kleborate_results.tsv > kleborate_mqc.tsv       ; fi
    if [ -f 'mlst_summary.tsv' ] ; then cut -f 1,3,4 mlst_summary.tsv > mlst_mqc.tsv            ; fi
    if [ -f 'plasmidfinder_result.tsv' ] ; then awk '{ print $1 "_" NR "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' plasmidfinder_result.tsv > plasmidfinder_mqc.tsv   ; fi
    if [ -f 'seqsero2_results.txt' ] ; then cut -f 1,4-10 seqsero2_results.txt > seqsero2_mqc.txt ; fi
    if [ -f 'serotypefinder_results.txt' ] ; then cut -f 1-6 serotypefinder_results.txt > serotypefinder_mqc.txt  ; fi

    multiqc !{params.multiqc_options} \
      --outdir multiqc \
      --cl_config "prokka_fn_snames: True"  \
      . \
      | tee -a $log_file
  '''
}

process core_genome_evaluation {
  tag         "Evaluating core genome"
  label       "process_single"
  publishDir  path: params.outdir, mode: 'copy', pattern: 'logs/*/*log'
  publishDir  path: params.outdir, mode: 'copy', pattern: 'core_genome_evaluation/core_genome_evaluation.csv'
  container   'quay.io/biocontainers/pandas:1.5.2'
  time        '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple file(fasta), file(summary), file(script)

  output:
  tuple file(fasta), env(num_samples), env(num_core_genes), emit: evaluation
  path "core_genome_evaluation/core_genome_evaluation.csv", emit: for_multiqc
  path "logs/${task.process}/*.log"                       , emit: log_files

  shell:
  def args = task.ext.args ?: ''
  """
    mkdir -p core_genome_evaluation logs/${task.process}
    log_file=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    python ${script} | tee -a \$log_file

    num_samples=\$(wc -l core_genome_evaluation.csv | awk '{print \$1}' )
    num_core_genes=\$(cut -f 3 core_genome_evaluation.csv -d "," | tail -n 1 | cut -f 1 -d "." )
    cp core_genome_evaluation.csv core_genome_evaluation/core_genome_evaluation.csv
  """
}

process flag {
  tag           "${meta.id}"
  label         "process_single"
  //publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.5.2'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), file(files)

  output:
  tuple val(meta), env(salmonella_flag)    , emit: salmonella_flag
  tuple val(meta), env(klebsiella_flag)    , emit: klebsiella_flag
  tuple val(meta), env(ecoli_flag)         , emit: ecoli_flag
  tuple val(meta), env(streppneu_flag)     , emit: streppneu_flag
  tuple val(meta), env(legionella_flag)    , emit: legionella_flag
  tuple val(meta), env(klebacin_flag)      , emit: klebacin_flag
  tuple val(meta), env(strepa_flag)        , emit: strepa_flag
  tuple val(meta), env(vibrio_flag)        , emit: vibrio_flag
  tuple val(meta), env(myco_flag)          , emit: myco_flag
  tuple val(meta), env(genus), env(species), emit: organism
  path "flag/*_flag.csv"                   , emit: collect
  path "logs/${task.process}/*.log"        , emit: log_files

  shell:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p flag logs/${task.process}
    log_file=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    if [ -f "${prefix}_fastani.csv" ]
    then 
      awk -F "," '{if (\$4 > 90) print \$0}' ${prefix}_fastani.csv > smaller_fastani.csv
      genus=\$(head   -n 2 ${prefix}_fastani.csv | tail -n 1 | cut -f 3 -d , | cut -f 1 -d "_")
      species=\$(head -n 2 ${prefix}_fastani.csv | tail -n 1 | cut -f 3 -d , | cut -f 2 -d "_")
    else
      touch smaller_fastani.csv
      genus="unknown"
      species="unknown"
    fi

    touch ${prefix}_reads_summary_kraken2.csv  ${prefix}.summary.mash.csv ${prefix}_blobtools.txt

    files="smaller_fastani.csv ${prefix}_reads_summary_kraken2.csv ${prefix}.summary.mash.csv ${prefix}_blobtools.txt"

    echo "Looking for Salmonella:" >> \$log_file
    salmonella_flag=''
    find_salmonella=\$(head -n 10 \$files | grep "Salmonella" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_salmonella" ] ; then salmonella_flag="found" ; else salmonella_flag="not" ; fi
    
    echo "Looking for E. coli and Shigella:" >> \$log_file
    ecoli_flag=''
    find_ecoli=\$(head -n 10 \$files | grep -e "Escherichia" -e "Shigella" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_ecoli" ] ; then ecoli_flag="found" ; else ecoli_flag="not" ; fi

    echo "Looking for Klebsiella:" >> \$log_file
    klebsiella_flag=''
    find_klebsiella=\$(head -n 10 \$files | grep -e "Klebsiella" -e "Enterobacter" -e "Serratia" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_klebsiella" ] ; then klebsiella_flag="found" ; else klebsiella_flag="not" ; fi
    
    echo "Looking for Strep A organisms:" >> \$log_file
    strepa_flag=''
    find_strepa=\$(head -n 10 \$files | grep "Streptococcus" | grep -e "pyogenes" -e "dysgalactiae" -e "anginosus" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_strepa" ] ; then strepa_flag='found' ; else strepa_flag='not' ; fi

    echo "Looking for Streptococcus pneumoniae organisms:" >> \$log_file
    streppneu_flag_flag=''
    find_streppneu=\$(head -n 10 \$files | grep "Streptococcus" | grep "pneumoniae" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_streppneu" ] ; then streppneu_flag='found' ; else streppneu_flag='not' ; fi
  
    echo "Looking for Legionella organisms:" >> \$log_file
    legionella_flag=''
    find_legionella=\$(head -n 10 \$files | grep "Legionella" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_legionella" ] ; then legionella_flag='found' ; else legionella_flag='not' ; fi

    echo "Looking for Vibrio organisms:" >> \$log_file
    vibrio_flag=''
    find_vibrio=\$(head -n 10 \$files | grep "Vibrio" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_vibrio" ] ; then vibrio_flag='found' ; else vibrio_flag='not' ; fi

    echo "Looking for Klebsiella or Acinetobacter:" >> \$log_file
    klebacin_flag=''
    if [ -n "\$find_klebsiella" ]
    then
      klebacin_flag='found'
    else
      find_acin=\$(head -n 10 \$files | grep "Acinetobacter" | tee -a \$log_file | head -n 1 )
      if [ -n "\$find_acin" ] ; then klebacin_flag='found' ; else klebacin_flag='not' ; fi
    fi

    echo "Looking for Mycobacterium/Mycobacteria"
    myco_flag=''
    find_myco=\$(head -n 10 \$files | grep "Mycobacteri" | tee -a \$log_file | head -n 1 )
    if [ -n "\$find_myco" ] ; then myco_flag='found' ; else myco_flag='not' ; fi

    if [ -z "\$genus" ]   ; then genus=unknown   ; fi
    if [ -z "\$species" ] ; then species=unknown ; fi

    echo "sample,genus,species,salmonella_flag,ecoli_flag,klebsiella_flag,klebacin_flag,myco_flag,strepa_flag,streppneu_flag,legionella_flag,vibrio_flag" > flag/${prefix}_flag.csv
    echo "${prefix},\$genus,\$species,\$salmonella_flag,\$ecoli_flag,\$klebsiella_flag,\$klebacin_flag,\$myco_flag,\$strepa_flag,\$streppneu_flag,\$legionella_flag,\$vibrio_flag" >> flag/${prefix}_flag.csv
  """
}

process json_convert {
  tag           "${prefix}"
  label         "process_single"
  // no publishDir
  container     'quay.io/biocontainers/pandas:1.5.2'
  maxForks      10
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), val(analysis), file(json), file(script)

  output:
  path "${analysis}/*_${analysis}_summary.csv", emit: collect

  shell:
  """
  mkdir -p ${analysis}

  python3 ${script} ${analysis} ${json} 
  """
}

process mash_err {
  tag           "${meta.id}"
  // no publishDir
  label         "process_single"
  container     'quay.io/biocontainers/pandas:1.5.2'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(error_file)

  output:
  path "summary/*_names.csv", emit: summary

  when:
  task.ext.when == null || task.ext.when

  shell:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """

  echo "whatever"

  exit 1
  """
}

process names {
  tag           "${meta.id}"
  // no publishDir
  label         "process_single"
  container     'quay.io/biocontainers/pandas:1.5.2'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  tuple val(meta), file(input)

  output:
  path "summary/*_names.csv", emit: collect

  when:
  task.ext.when == null || task.ext.when

  shell:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p summary

  echo "sample,file,version" > summary/${prefix}_names.csv
  echo "${prefix},${input},${workflow.manifest.version}" >> summary/${prefix}_names.csv
  """
}

process references {
  tag       "Preparing references"
  // no publishDir
  label     "process_single"
  container 'quay.io/uphl/grandeur_ref:20240124'
  time      '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  output:
  path "ref/*", emit: fastas

  when:
  task.ext.when == null || task.ext.when

  shell:
  """
  mkdir ref

  cp /ref/* ref/.
  """
}

process species {
  tag           "Creating list of species"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.5.2'
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  
  input:
  file(results)

  output:
  path "datasets/species_list.txt", emit: species                                    

  when:
  task.ext.when == null || task.ext.when

  shell:
  """
    mkdir -p datasets

    if [ -f "blobtools_species.txt" ]
    then 
      cut -f 2 blobtools_species.txt >> species.txt
    fi

    if [ -f "kraken2_summary.csv" ]
    then 
      cut -f 7 -d , kraken2_summary.csv >> species.txt
    fi

    if [ -f "mash_summary.csv" ]
    then
      cut -f 7 -d , mash_summary.csv | tail -n+2 >> species.txt
    fi

    grep -v no-hit species.txt | grep -v undef | grep -v name | grep "_" | sort | uniq > datasets/species_list.txt
  """
}

process summary {
  tag           "Creating summary files"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/pandas:1.5.2'
  label         "process_single"
  time          '10m'
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(input)

  output:
  path "grandeur_summary.tsv"                 , emit: summary_tsv
  path "grandeur_summary.txt"                 , emit: summary_txt
  path "summary/grandeur_extended_summary.tsv", emit: extended_tsv
  path "summary/grandeur_extended_summary.txt", emit: extended_txt

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args = task.ext.args ?: ''
  """
    mkdir -p summary

    python summary.py
  """
}

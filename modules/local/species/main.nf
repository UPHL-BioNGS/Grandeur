process SPECIES {
  tag           "Creating list of species"
  label         "process_single"
  container     'staphb/pandas:3.0.1'
  
  input:
  path results // Changed from file(results) for modern Nextflow syntax

  output:
  path "datasets/species_list.txt", emit: species, optional: true
  path "datasets/accessions_list.txt", emit: accessions, optional: true                                 

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    mkdir -p datasets
    
    # Ensure files exist so grep doesn't crash if no summaries are found
    touch species.txt accessions.txt

    # kraken2
    if [ -f "kraken2_summary.csv" ]; then 
      cut -f 7 -d , kraken2_summary.csv >> species.txt
    fi

    # mash dist
    if [ -f "mashdist_summary.csv" ]; then
      cut -f 7 -d , mashdist_summary.csv | tail -n+2 >> species.txt
    fi

    # mash screen
    if [ -f "mashscreen_summary.csv" ]; then
      cut -f 7 -d , mashscreen_summary.csv | tail -n+2 >> species.txt
    fi

    # spestimator
    if [ -f "spestimator_summary.tsv" ]; then
      cut -f 3 -d , spestimator_summary.tsv | awk '{print \$1 "_" \$2}' | tail -n+2 >> species.txt
      cut -f 4 -d , spestimator_summary.tsv | tail -n+2 >> accessions.txt
    fi  

    # sylph
    if [ -f "sylph_summary.csv" ]; then
      cat sylph_summary.csv >> accessions.txt
    fi
    
    grep "_" species.txt | sed 's/^_//g' | grep -v "_sp\\." |  grep -v "_sp\$" | sort | uniq > datasets/species_list.txt
    grep "_" accessions.txt | sort | uniq > datasets/accessions_list.txt
  """
}
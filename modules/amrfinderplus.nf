process amrfinderplus {
  tag "${sample}"
  label "medcpus"

  when:
  params.contig_processes =~ /amrfinderplus/

  input:
  tuple val(sample), file(contigs), val(genus), val(species)

  output:
  path "ncbi-AMRFinderplus/${sample}_amrfinder_plus.txt"                , emit: collect
  tuple val(sample), env(amr_genes)                                     , emit: amr_genes
  tuple val(sample), env(virulence_genes)                               , emit: vir_genes
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p ncbi-AMRFinderplus logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    organism_options=(Acinetobacter_baumannii
    Campylobacter
    Enterococcus_faecalis
    Enterococcus_faecium
    Escherichia:Shigella
    Klebsiella
    Salmonella
    Staphylococcus_aureus
    Staphylococcus_pseudintermedius
    Streptococcus_agalactiae
    Streptococcus_pneumoniae
    Streptococcus_pyogenes
    Vibrio_cholerae)

    organism=$(history -p ${organism_options[@]} | grep -i !{genus} | grep -i !{species} | head -n 1 )
    if [ -z "$organism" ] ; then organism=$(history -p ${organism_options[@]} | grep -i !{genus} | head -n 1 | cut -f 1 -d ":" ) ; fi
    if [ -n "$organism" ]
    then
      organism_check="--organism $organism"
      echo "Mash result of !{genus} !{species} matched with $organism" >> $log_file
    else
      organism_check=''
      echo "Mash result of !{genus} !{species} did not match any of the organisms" >> $log_file
    fi

    amrfinder !{params.amrfinderplus_options} \
      --nucleotide !{contigs} \
      --threads !{task.cpus} \
      --name !{sample} \
      --output ncbi-AMRFinderplus/!{sample}_amrfinder_plus.txt \
      $organism_check \
      --plus \
      2>> $err_file >> $log_file

    amr_genes=$(cut -f 7 ncbi-AMRFinderplus/!{sample}_amrfinder_plus.txt | tail +2 | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    virulence_genes=$(grep "VIRULENCE" ncbi-AMRFinderplus/!{sample}_amrfinder_plus.txt | cut -f 7 | sort | uniq | tr '\\n' ',' | sed 's/,$//g' )
    if [ -z "$amr_genes" ] ; then amr_genes="none" ; fi
    if [ -z "$virulence_genes" ] ; then virulence_genes="none" ; fi
  '''
}

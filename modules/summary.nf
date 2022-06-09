process summary {
  tag "${sample}"

  when:
  params.fastq_processes =~ /summary/ || params.contig_processes =~ /summary/

  input:
  tuple val(sample), file(file),
    // de_novo_alignment
    val(bbduk_phix_results),
    val(fastp_results),

    // fastq_information
    val(fastqc_1_results),
    val(fastqc_2_results),
    val(cg_avrl_results),
    val(cg_quality_results),
    val(cg_cov_results),
    val(ref_genome_length),
    val(shigatyper_predictions),
    val(shigatyper_cadA),
    val(kraken2_top_hit),
    val(kraken2_top_perc),
    val(kraken2_top_reads),

    // mash
    val(mash_genome_size_results),
    val(mash_coverage_results),
    val(mash_genus_results),
    val(mash_species_results),
    val(mash_full_results),
    val(mash_pvalue_results),
    val(mash_distance_results),

    // contig_information
    val(seqsero2_profile_results),
    val(seqsero2_serotype_results),
    val(seqsero2_contamination_results),
    val(serotypefinder_results_o),
    val(serotypefinder_results_h),
    val(kraken2_top_hit_contigs),
    val(kraken2_top_perc_contigs),
    val(kraken2_top_reads_contigs),
    val(plasmidfinder_hits),
    val(quast_gc_results),
    val(quast_contigs_results),
    val(quast_N50_contigs_results),
    val(quast_length_results),
    val(kleborate_score),
    val(kleborate_mlst),
    val(amr_genes),
    val(virulence_genes),
    val(fastani_ref_results),
    val(fastani_ani_results),
    val(fastani_fragment),
    val(fastani_total),
    val(mlst_results),

    // blobtools
    val(blobtools_species_results),
    val(blobtools_perc_results)

  output:
  path "summary/${sample}.summary.txt"                                  , emit: summary_files_txt
  path "summary/${sample}.summary.tsv"                                  , emit: summary_files_tsv
  path "logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}" , emit: log

  shell:
  '''
    mkdir -p summary logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "container : !{task.container}" >> $log_file
    echo "Nextflow command : " >> $log_file
    cat .command.sh >> $log_file

    sample_id_split=($(echo !{sample} | sed 's/-/ /g' | sed 's/_/ /g' ))
    if [ "${#sample_id_split[@]}" -ge "5" ]
    then
      sample_id="${sample_id_split[0]}-${sample_id_split[1]}"
    elif [ "${#sample_id_split[@]}" -eq "4" ]
    then
      sample_id=${sample_id_split[0]}
    else
      sample_id=!{sample}
    fi

    header="sample_id;sample"
    result="$sample_id;!{sample}"

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep fastp)" ]
    then
      header="$header;reads_passed"
      result="$result;!{fastp_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep bbduk)" ]
    then
      header="$header;phix_reads_removed"
      result="$result;!{bbduk_phix_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep fastqc)" ]
    then
      header="$header;fastqc_1_reads;fastqc_2_reads"
      result="$result;!{fastqc_1_results};!{fastqc_2_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep mash)" ]
    then
      header="$header;mash_genome_size;mash_coverage;mash_genus;mash_species;mash_full;mash_pvalue;mash_distance"
      result="$result;!{mash_genome_size_results};!{mash_coverage_results};!{mash_genus_results};!{mash_species_results};!{mash_full_results};!{mash_pvalue_results};!{mash_distance_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep fastani)" ]
    then
      header="$header;fastani_ref_top_hit;fastani_ani_score;fastani_fragment_mappings;fastani_total"
      result="$result;!{fastani_ref_results};!{fastani_ani_results};!{fastani_fragment};!{fastani_total}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep quast)" ]
    then
      header="$header;quast_gc_%;quast_contigs;quast_N50;quast_length"
      result="$result;!{quast_gc_results};!{quast_contigs_results};!{quast_N50_contigs_results};!{quast_length_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep cg_pipeline)" ]
    then
      header="$header;cg_average_read_length;cg_average_quality;cg_coverage;ref_genome_length"
      result="$result;!{cg_avrl_results};!{cg_quality_results};!{cg_cov_results};!{ref_genome_length}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep seqsero2)" ]
    then
      header="$header;seqsero2_profile;seqsero2_serotype;seqsero2_contamination"
      result="$result;!{seqsero2_profile_results};!{seqsero2_serotype_results};!{seqsero2_contamination_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep serotypefinder)" ]
    then
      header="$header;serotypefinder_o_group;serotypefinder_h_group"
      result="$result;!{serotypefinder_results_o};!{serotypefinder_results_h}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep kleborate)" ]
    then
      header="$header;kleborate_score;kleborate_mlst"
      result="$result;!{kleborate_score};!{kleborate_mlst}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep amrfinderplus)" ]
    then
      header="$header;amr_genes;virulence_genes"
      result="$result;!{amr_genes};!{virulence_genes}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep blobtools)" ] && [ "!{params.blast_db}" != "false" ]
    then
      header="$header;blobtools_top_species;blobtools_percentage"
      result="$result;!{blobtools_species_results};!{blobtools_perc_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep kraken2)" ] && [ "!{params.kraken2_db}" != "false" ]
    then
      header="$header;kraken2_top_species_reads;kraken2_num_reads;kraken2_percentage_reads;kraken2_top_species_contigs;kraken2_num_contigs;kraken2_percentage_contigs"
      result="$result;!{kraken2_top_hit};!{kraken2_top_reads};!{kraken2_top_perc};!{kraken2_top_hit_contigs};!{kraken2_top_reads_contigs};!{kraken2_top_perc_contigs}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep mlst)" ]
    then
      header="$header;mlst"
      result="$result;!{mlst_results}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep shigatyper)" ]
    then
      header="$header;shigatyper_predictions;shigatyper_cadA"
      result="$result;!{shigatyper_predictions};!{shigatyper_cadA}"
    fi

    if [ -n "$(echo !{params.fastq_processes} !{params.contig_processes} !{params.phylogenetic_processes} | grep plasmidfinder)" ]
    then
      header="$header;plasmidfinder_hits"
      result="$result;!{plasmidfinder_hits}"
    fi

    echo $header > summary/!{sample}.summary.txt
    echo $result >> summary/!{sample}.summary.txt

    cat summary/!{sample}.summary.txt | tr '\t' '_' | tr ';' '\t' > summary/!{sample}.summary.tsv
  '''
}

#!/usr/bin/env nextflow

println("Currently using the Peaks workflow for use with annotated contig files\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210830")
println("")

// TODO : add ete3 or some tree building software

params.outdir = workflow.launchDir + '/peaks'
println("The files and directory for results is " + params.outdir)

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")

params.contigs = workflow.launchDir + '/contigs'
params.from_contigs = true
contigs = params.from_contigs
  ? Channel
      .fromPath("${params.contigs}/*.{fa,fasta,fna}", type: 'file')
      .map { file -> tuple(file.baseName, file) }
      .view { "fasta file : $it" }
  : Channel.empty()

params.gff = workflow.launchDir + '/gff'
params.from_gff = true
local_gffs = params.from_gff
  ? Channel.fromPath("${params.gff}/*.gff", type: 'file').view { "gff file : $it" }
  : Channel.empty()

params.kraken2 = false
params.kraken2_db = '/kraken2_db'
local_kraken2 = params.kraken2
  ? Channel
    .fromPath(params.kraken2_db, type:'dir')
    .view { "Local kraken2 database : $it" }
    .ifEmpty{
      println("No kraken2 database was found at ${params.kraken2_db}")
      println("Set 'params.kraken2_db' to directory with kraken2 database")
      exit 1
    }
  : Channel.empty()

params.prokka = true
params.prokka_options = ''
params.center = 'STAPHB'
params.mincontiglen = 500
process prokka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus
  container 'staphb/prokka:latest'

  when:
  params.prokka

  input:
  tuple val(sample), file(contigs) from contigs

  output:
  file("prokka/${sample}/${sample}.{err,faa,ffn,fna,fsa,gbk,log,sqn,tbl,tsv,txt}")
  file("prokka/${sample}/${sample}.gff") into prokka_gffs
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    prokka -v >> $log_file

    prokka !{params.prokka_options} \
      --cpu !{task.cpus} \
      --compliant \
      --centre !{params.center} \
      --mincontiglen !{params.mincontiglen} \
      --outdir prokka/!{sample} \
      --locustag locus_tag \
      --prefix !{sample} \
      --force !{contigs} \
      2>> $err_file >> $log_file
  '''
}

prokka_gffs
  .concat(local_gffs)
  .ifEmpty{
    println("No gff files were found at ${params.gff}")
    println("No contig or fasta files ending with '.fa', '.fna', or '.fasta' were found at ${params.contigs}")
    println("Set 'params.gff' to directory with gff files")
    println("Set 'params.contigs' to directory with fastas")
    exit 1
  }
  .set { gffs }

params.roary = true
params.roary_options = ''
if (params.kraken2) {
  process roary_qc {
    publishDir "${params.outdir}", mode: 'copy'
    echo false
    cpus params.maxcpus
    container 'staphb/roary:latest'

    when:
    params.roary

    input:
    file(contigs) from gffs.collect()
    path(local_kraken2) from local_kraken2

    output:
    file("roary/*")
    file("roary/fixed_input_files/*")
    file("roary/core_gene_alignment.aln") into roary_core_genome_iqtree, roary_core_genome_snp_dists
    file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p logs/!{task.process}
      log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      roary -a >> $log_file

      roary !{params.roary_options} \
        -p !{task.cpus} \
        -f roary \
        -e -n \
        -qc -k !{local_kraken2} \
        *.gff \
        2>> $err_file >> $log_file
    '''
  }
} else {
  process roary {
    publishDir "${params.outdir}", mode: 'copy'
    echo false
    cpus params.maxcpus
    container 'staphb/roary:latest'

    when:
    params.roary

    input:
    file(contigs) from gffs.collect()

    output:
    file("${task.process}/*")
    file("${task.process}/fixed_input_files/*")
    file("${task.process}/core_gene_alignment.aln") into roary_core_genome_iqtree, roary_core_genome_snp_dists
    file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p logs/!{task.process}
      log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      roary -a >> $log_file

      roary !{params.roary_options} \
        -p !{task.cpus} \
        -f !{task.process} \
        -e -n \
        *.gff \
        2>> $err_file >> $log_file
    '''
  }
}

params.iqtree = true
params.iqtree_options = '-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000'
process iqtree {
  publishDir "${params.outdir}", mode: 'copy'
  echo false
  cpus params.maxcpus
  container 'staphb/iqtree:latest'

  when:
  params.iqtree

  input:
  file(msa) from roary_core_genome_iqtree

  output:
  file("${task.process}/iqtree{.ckp.gz,.treefile,.iqtree,.log,.splits.nex}")
  file("${task.process}/iqtree.contree") into treefile
  file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    iqtree -v >> $log_file

    iqtree !{params.iqtree_options} \
      -s !{msa} \
      -pre !{task.process}/!{task.process} \
      -nt AUTO \
      -ntmax !{task.cpus} \
      2>> $err_file >> $log_file
  '''
}

// WARNING : THIS CONTAINER DOESN'T EXIST
params.ete3 = false
params.ete3_options = ''
process ete3 {
  publishDir "${params.outdir}", mode: 'copy'
  cpus 1
  //container 'staphb/ete3:latest'
  //container 'docker://quay.io/biocontainers/ete3:3.1.2'

  when:
  params.ete3

  input:
  file(newick) from treefile

  output:
  file("${task.process}/tree.svg")
  file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "ETE3 version : $(ete3 version | head -n 1 )" | tee -a $log_file

    ete3 view --image !{task.process}/tree.svg -t !{newick}
  '''
}

params.snp_dists = true
params.snp_dists_options = ''
process snp_dists {
  publishDir "${params.outdir}", mode: 'copy'
  echo false
  cpus 1
  container 'staphb/snp-dists:latest'

  when:
  params.snp_dists

  input:
  file(contigs) from roary_core_genome_snp_dists

  output:
  file("${task.process}/snp_matrix.txt")
  file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    snp-dists -v >> $log_file

    snp-dists !{params.snp_dists_options} \
      !{contigs} \
      2>> $err_file \
      > !{task.process}/snp_matrix.txt
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}

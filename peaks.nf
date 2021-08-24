#!/usr/bin/env nextflow

println("Currently using the Peaks workflow for use with annotated contig files\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210830")
println("")

// TODO add ETE3

params.outdir = workflow.launchDir + '/peaks'
println("The files and directory for results is " + params.outdir)

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

params.contigs = workflow.launchDir + '/contigs'
params.from_contigs = true
contigs = params.from_contigs
  ? Channel
      .fromPath("${params.contigs}/*.{fa,fasta,fna}", type: 'file')
      .map { file -> tuple(file.baseName, file) }
      .view { "Contigs or Fasta files : $it" }
  : Channel.empty()

params.gff = workflow.launchDir + '/gff'
params.from_gff = true
local_gffs = params.from_gff
  ? Channel.fromPath("${params.gff}/*.gff", type: 'file').view { "gff files : $it" }
  : Channel.empty()

params.kraken2 = false
params.local_kraken2 = '/kraken2_db'
local_kraken2 = params.kraken2
  ? Channel
    .fromPath(params.local_kraken2, type:'dir')
    .view { "Local kraken2 database : $it" }
    .ifEmpty{
      println("No kraken2 database was found at ${params.local_kraken2}")
      println("Set 'params.local_kraken2' to directory with kraken2 database")
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
  process roary {
    publishDir "${params.outdir}", mode: 'copy'
    echo false
    cpus params.maxcpus
    container 'staphb/roary:latest'

    when:
    params.roary

    input:
    file(contigs), dir(kraken2) from gffs.collect()
    dir(local_kraken2) from local_kraken2

    output:
    file("output") into roary_core_genome_iqtree, roary_core_genome_snp_dists
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      roary -a >> $log_file

      roary !{params.roary_options} \
        -p !{task.cpus} \
        -f !{task.process} \
        -e -n \
        -qc -k !{kraken2} \
        *.gff \
        2>> $err_file >> $log_file

      exit 1
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
    file(contigs), dir(kraken2) from gffs.collect()

    output:
    file("output") into roary_core_genome_iqtree, roary_core_genome_snp_dists
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
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

      exit 1
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
  file("output")
  file("output") into treefile
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

    exit 1
  '''
}

params.ete3 = true
params.ete3_options = ''
process ete3 {
  publishDir "${params.outdir}", mode: 'copy'
  echo false
  cpus params.maxcpus
  container 'reslp/ete3:3.1.1'
  // barikan/ete3_with_qt:3.1.1
  // bhntools/dockerete:latest

  when:
  params.ete3

  input:
  file(treefile) from treefile

  output:
  file("output")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    ete3 -v >> $log_file

    ete3

    exit 1
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
  file("${task.process}/results.txt")
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
      > !{task.process}/results.txt

    exit 1
  '''
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}

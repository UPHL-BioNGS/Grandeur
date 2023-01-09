include { amrfinderplus }  from '../modules/amrfinderplus'  addParams(params)
include { fastqc }         from '../modules/fastqc'         addParams(params)
include { fastqscan }      from '../modules/fastqscan'      addParams(params)
include { flag }           from '../modules/grandeur'       addParams(params)
include { kleborate }      from '../modules/kleborate'      addParams(params)
include { mlst }           from '../modules/mlst'           addParams(params)
include { plasmidfinder }  from '../modules/plasmidfinder'  addParams(params)
include { quast }          from '../modules/quast'          addParams(params)
include { seqsero2 }       from '../modules/seqsero2'       addParams(params)
include { serotypefinder } from '../modules/serotypefinder' addParams(params)
include { shigatyper }     from '../modules/shigatyper'     addParams(params)
include { size }           from '../modules/grandeur'       addParams(params)

workflow information {
  take:
    ch_reads
    ch_contigs
    ch_species
    ch_size

  main:
    fastqc(ch_reads)
    
    mlst(ch_contigs)
    quast(ch_contigs)
    plasmidfinder(ch_contigs)

    flag(ch_species.groupTuple(by: 0))

    size(ch_size)
    fastqscan(ch_reads.join(size.out.size, by: 0))
    
    kleborate(ch_contigs.join(flag.out.klebsiella_flag, by:0))
    seqsero2(ch_contigs.join(flag.out.salmonella_flag,  by:0))
    serotypefinder(ch_contigs.join(flag.out.ecoli_flag, by:0))
    shigatyper(ch_contigs.join(flag.out.ecoli_flag,     by:0))
    
    amrfinderplus(ch_contigs.join(flag.out.organism,    by:0))
    
    fastqc.out.collect
      .collectFile(name: "fastqc_summary.csv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/fastqc")
      .set{ fastqc_summary }

    fastqscan.out.collect
      .collectFile(name: "fastqscan_summary.csv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/fastq-scan")
      .set{ fastqscan_summary }

    mlst.out.collect
      .collectFile(name: "mlst_summary.csv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/mlst")
      .set{ mlst_summary }

    quast.out.collect
      .collectFile(name: "quast_report.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/quast")
      .set{ quast_summary }

    plasmidfinder.out.collect
      .collectFile(name: "plasmidfinder_result.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/plasmidfinder")
      .set{ plasmidfinder_summary }

    kleborate.out.collect
      .collectFile(name: "kleborate_results.txt",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kleborate")
      .set{ kleborate_summary }

    seqsero2.out.collect
      .collectFile(name: "seqsero2_results.txt",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/seqsero2")
      .set{ seqsero2_summary }

    serotypefinder.out.collect
      .collectFile(name: "serotypefinder_results.txt",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/serotypefinder")
      .set{ serotypefinder_summary }

    shigatyper.out.collect
      .collectFile(name: "shigatyper_results.txt",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigatyper")
      .set{ shigatyper_summary }

    amrfinderplus.out.collect
      .collectFile(name: "amrfinderplus.txt",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/ncbi-AMRFinderplus")
      .set{ amrfinderplus_summary }

    amrfinderplus_summary
      .mix(fastqc_summary)
      .mix(fastqscan_summary)
      .mix(kleborate_summary)
      .mix(mlst_summary)
      .mix(plasmidfinder_summary)
      .mix(quast_summary)
      .mix(seqsero2_summary)
      .mix(serotypefinder_summary)
      .mix(shigatyper_summary)
      .set { for_summary }

    fastqc.out.for_multiqc
      .mix(quast.out.for_multiqc)
      .set { for_multiqc }

  emit:
    for_summary = for_summary.collect()
    for_multiqc = for_multiqc.collect()
}

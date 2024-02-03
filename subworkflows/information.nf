include { amrfinderplus }  from '../modules/local/amrfinderplus'  addParams(params)
include { drprg }          from '../modules/local/drprg'          addParams(params)
include { elgato }         from '../modules/local/elgato'         addParams(params)
include { emmtyper }       from '../modules/local/emmtyper'       addParams(params)
include { fastqc }         from '../modules/local/fastqc'         addParams(params)
include { flag }           from '../modules/local/local'          addParams(params)
include { json_convert }   from '../modules/local/local'          addParams(params)
include { kaptive }        from '../modules/local/kaptive'        addParams(params)
include { kleborate }      from '../modules/local/kleborate'      addParams(params)
include { legsta }         from '../modules/local/legsta'         addParams(params)
include { mlst }           from '../modules/local/mlst'           addParams(params)
include { mykrobe }        from '../modules/local/mykrobe'        addParams(params)
include { pbptyper }       from '../modules/local/pbptyper'       addParams(params)
include { plasmidfinder }  from '../modules/local/plasmidfinder'  addParams(params)
include { quast }          from '../modules/local/quast'          addParams(params)
include { seqsero2 }       from '../modules/local/seqsero2'       addParams(params)
include { serotypefinder } from '../modules/local/serotypefinder' addParams(params)
include { shigatyper }     from '../modules/local/shigatyper'     addParams(params)

workflow information {
  take:
    ch_reads
    ch_contigs
    ch_flag
    summfle_script
    jsoncon_script

  main:
    for_multiqc = Channel.empty()
    ch_versions = Channel.empty()

    // fastq files only
    if ( params.sample_sheet || params.reads ) {
      fastqc(ch_reads)

      for_multiqc = for_multiqc.mix(fastqc.out.for_multiqc)

      fastqc.out.collect
        .collectFile(name: "fastqc_summary.csv",
          keepHeader: true,
          sort: { file -> file.text },
          storeDir: "${params.outdir}/fastqc")
        .set{ fastqc_summary }

      ch_versions = ch_versions.mix(fastqc.out.versions)

    } else {
      fastqc_summary = Channel.empty()
    }

    // contigs
    mlst(ch_contigs.combine(summfle_script))
    quast(ch_contigs)
    plasmidfinder(ch_contigs.combine(summfle_script))

    // species specific
    flag(ch_flag.groupTuple())

    amrfinderplus(ch_contigs.join(flag.out.organism,    by:0))
    drprg(ch_contigs.join(flag.out.myco_flag,           by:0))
    emmtyper(ch_contigs.join(flag.out.strepa_flag,      by:0).combine(summfle_script)) 
    //kaptive(ch_contigs.join(flag.out.klebacin_flag,     by:0))      
    kleborate(ch_contigs.join(flag.out.klebsiella_flag, by:0).combine(summfle_script))
    elgato(ch_contigs.join(flag.out.legionella_flag,    by:0))
    mykrobe(ch_contigs.join(flag.out.myco_flag,         by:0))
    pbptyper(ch_contigs.join(flag.out.streppneu_flag,   by:0))
    seqsero2(ch_contigs.join(flag.out.salmonella_flag,  by:0))
    serotypefinder(ch_contigs.join(flag.out.ecoli_flag, by:0).combine(summfle_script))
    shigatyper(ch_contigs.join(flag.out.ecoli_flag,     by:0).combine(summfle_script))

    json_convert(drprg.out.json.combine(jsoncon_script))

    amrfinderplus.out.collect
      .collectFile(name: "amrfinderplus.txt",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/ncbi-AMRFinderplus")
      .set{ amrfinderplus_summary }

    json_convert.out.collect
      .filter( /drprg/ )
      .collectFile(name: "drprg_summary.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/drprg")
      .set{ drprg_summary }

    elgato.out.collect
      .collectFile(name: "elgato_summary.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/elgato")
      .set{ elgato_summary }

    emmtyper.out.collect
      .collectFile(name: "emmtyper_summary.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/emmtyper")
      .set{ emmtyper_summary }

    flag.out.collect
      .collectFile(name: "flag_summary.csv",
        keepHeader: true,
        sort: {file -> file.text },
        storeDir: "${params.outdir}/flag")
      .set { flag_summary }

      // kaptive.out.collect
      //   .collectFile(name: "kaptive_summary.csv",
      //     keepHeader: true,
      //     sort: { file -> file.text },
      //     storeDir: "${params.outdir}/kaptive")
      //   .set{ kaptive_summary }

    kleborate.out.collect
      .collectFile(name: "kleborate_results.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kleborate")
      .set{ kleborate_summary }

    mlst.out.collect
      .collectFile(name: "mlst_summary.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/mlst")
      .set{ mlst_summary }

    mykrobe.out.collect
      .collectFile(name: "mykrobe_summary.csv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/mykrobe")
      .set{ mykrobe_summary }

    quast.out.collect
      .collectFile(name: "quast_report.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/quast")
      .set{ quast_summary }

    pbptyper.out.collect
      .collectFile(name: "pbptyper_summary.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/pbptyper")
      .set{ pbptyper_summary }

    plasmidfinder.out.collect
      .collectFile(name: "plasmidfinder_result.tsv",
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/plasmidfinder")
      .set{ plasmidfinder_summary }

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

    // TODO:
    // check all are in mutliqc
    // check all are in summary file
    // check kleborate getting amr
    // known issues : emmtyper, kleborate, el gato, pbp, drprg


    amrfinderplus_summary
      .mix(drprg_summary)
      .mix(elgato_summary)
      .mix(emmtyper_summary)
      .mix(fastqc_summary)
      //.mix(kaptive_summary)
      .mix(kleborate_summary)
      .mix(mlst_summary)
      .mix(mykrobe_summary)
      .mix(pbptyper_summary)
      .mix(plasmidfinder_summary)
      .mix(quast_summary)
      .mix(seqsero2_summary)
      .mix(serotypefinder_summary)
      .mix(shigatyper_summary)
      .set { for_summary }

    ch_versions
      .mix(amrfinderplus.out.versions)
      .mix(drprg.out.versions)
      .mix(elgato.out.versions)
      .mix(emmtyper.out.versions)
      //.mix(kaptive.out.versions)
      .mix(kleborate.out.versions)
      .mix(mlst.out.versions)
      .mix(mykrobe.out.versions)
      .mix(pbptyper.out.versions)
      .mix(plasmidfinder.out.versions)
      .mix(quast.out.versions)
      .mix(seqsero2.out.versions)
      .mix(serotypefinder.out.versions)
      .mix(shigatyper.out.versions)
      .set { for_versions }

  emit:
    for_summary = for_summary.collect()
    for_multiqc = for_multiqc.mix(quast.out.for_multiqc).collect()
    versions    = for_versions
}


include { DRPRG }          from '../../../modules/local/drprg'
include { ELGATO }         from '../../../modules/local/elgato'
include { EMMTYPER }       from '../../../modules/local/emmtyper'
include { JSON_CONVERT }   from '../../../modules/local/json_convert'
include { KAPTIVE }        from '../../../modules/local/kaptive'
include { KLEBORATE }      from '../../../modules/local/kleborate'
include { MENINGOTYPE }    from '../../../modules/local/meningotype'
include { MYKROBE }        from '../../../modules/local/mykrobe'
include { NGMASTER }       from '../../../modules/local/ngmaster'
include { PBPTYPER }       from '../../../modules/local/pbptyper'
include { SEQSERO2 }       from '../../../modules/local/seqsero2'
include { SEROTYPEFINDER } from '../../../modules/local/serotypefinder'
include { SHIGAPASS }      from '../../../modules/local/shigapass'


workflow SUBTYPING {
  take:
    ch_myco
    ch_gas
    ch_kleb
    ch_legionella
    ch_strep
    ch_salmonella
    ch_ecoli
    ch_vibrio
    ch_gc
    summfle_script
    jsoncon_script

  main:
    ch_summary  = Channel.empty()
    ch_versions = Channel.empty()

    log.info "Running subtyping analysis. This workflow will perform in silico subtyping of assemblies with a variety of tools, depending on the species of interest."
    log.info "Current species-specific subtyping tools are for Mycobacterium, Streptococcus pneumoniae, Salmonella, Escherichia coli, Vibrio, Neisseria meningitidis, Neisseria gonorrhoeae, Klebsiella, Legionella, and Streptococcus pyogenes."
    log.info "If a desired sub-typing tool is not included here, please contact the developers or submit an issue on GitHub at https://github.com/UPHL-BioNGS/Grandeur/issues"

    log.info "DR_PRG is a tool for predicting the drug resistance phenotype of Mycobacterium tuberculosis from whole genome sequencing data."
    DRPRG(ch_myco.filter{it})

    JSON_CONVERT(DRPRG.out.json.combine(jsoncon_script))

    JSON_CONVERT.out.collect
      .filter( ~/.*drprg.tsv/ )
      .collectFile(name: 'drprg_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/drprg")
      .set{ drprg_summary }

    ch_summary  = ch_summary.mix(drprg_summary)
    ch_versions = ch_versions.mix(DRPRG.out.versions.first())

    log.info "EMMtyper is a tool for in silico emm typing of Streptococcus pyogenes assemblies."
    EMMTYPER(ch_gas.filter{it}.combine(summfle_script)) 

    EMMTYPER.out.collect
      .collectFile(name: 'emmtyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/emmtyper")
      .set{ emmtyper_summary }

    ch_summary  = ch_summary.mix(emmtyper_summary)
    ch_versions = ch_versions.mix(EMMTYPER.out.versions.first())

    log.info "KAPTIVE is used for in silico K and O locus typing of Vibrio assemblies, with a focus on Vibrio parahaemolyticus."
    KAPTIVE(ch_vibrio.filter{it})      

    KAPTIVE.out.collect
      .collectFile(name: 'kaptive_summary.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kaptive")
      .set{ kaptive_summary }
    
    ch_summary  = ch_summary.mix(kaptive_summary)
    ch_versions = ch_versions.mix(KAPTIVE.out.versions.first())

    log.info "Kleborate is a tool for in silico subtyping of Klebsiella assemblies, including species assignment, multi-locus sequence typing, K and O locus typing, and detection of virulence and AMR genes."
    KLEBORATE(ch_kleb.filter{it}.combine(summfle_script))

    KLEBORATE.out.collect
      .collectFile(name: 'kleborate_results.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kleborate")
      .set{ kleborate_summary }
    
    ch_summary  = ch_summary.mix(kleborate_summary)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())

    log.info "EL GATO is a tool for in silico subtyping of Legionella assemblies, including species assignment, multi-locus sequence typing, and detection of virulence genes."
    ELGATO(ch_legionella.filter{it})

    ELGATO.out.collect
      .collectFile(name: 'elgato_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/elgato")
      .set{ elgato_summary }

    ch_summary = ch_summary.mix(elgato_summary)
    ch_versions = ch_versions.mix(ELGATO.out.versions.first())

    log.info "MYKROBE is a tool for predicting the drug resistance phenotype of Mycobacterium from whole genome sequencing data."
    MYKROBE(ch_myco.filter{it})

    MYKROBE.out.collect
      .collectFile(name: 'mykrobe_summary.csv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/mykrobe")
      .set{ mykrobe_summary }

    ch_summary  = ch_summary.mix(mykrobe_summary)
    ch_versions = ch_versions.mix(MYKROBE.out.versions.first())

    log.info "MENINGOTYPE is a tool for in silico subtyping of Neisseria meningitidis assemblies, including multi-locus sequence typing and PorA and FetA subtyping."
    MENINGOTYPE(ch_gc.filter{it})

    MENINGOTYPE.out.files
      .collectFile(name: 'meningotype_summary.tsv',
        keepHeader: true,
        sort: {file -> file.text },
        storeDir: "${params.outdir}/meningotype")
      .set{ meningotype_summary }

    ch_summary  = ch_summary.mix(meningotype_summary)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions.first())

    log.info "NGMASTER is a tool for in silico NG-MAST typing of Neisseria gonorrhoeae assemblies."
    NGMASTER(ch_gc.filter{it}.combine(summfle_script))

    NGMASTER.out.collect
      .collectFile(name: 'ngmaster_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/ngmaster")
      .set{ ngmaster_summary }

    ch_summary  = ch_summary.mix(ngmaster_summary)
    ch_versions = ch_versions.mix(NGMASTER.out.versions.first())

    log.info "PBPTYPER is a tool for in silico penicillin binding protein (PBP) typer of Streptococcus pneumoniae assemblies."
    PBPTYPER(ch_strep.filter{it})

    PBPTYPER.out.collect
      .collectFile(name: 'pbptyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/pbptyper")
      .set{ pbptyper_summary }

    ch_summary  = ch_summary.mix(pbptyper_summary)
    ch_versions = ch_versions.mix(PBPTYPER.out.versions.first())

    log.info "SEQSERO2 is a tool for in silico serotyping of Salmonella assemblies."
    SEQSERO2(ch_salmonella.filter{it})

    SEQSERO2.out.collect
      .collectFile(name: 'seqsero2_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/seqsero2")
      .set{ seqsero2_summary }

    ch_summary  = ch_summary.mix(seqsero2_summary)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions.first())

    log.info "SEROTYPERFINDER is a tool for in silico serotyping of Escherichia coli assemblies."
    SEROTYPEFINDER(ch_ecoli.filter{it}.combine(summfle_script))

    SEROTYPEFINDER.out.collect
      .collectFile(name: 'serotypefinder_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/serotypefinder")
      .set{ serotypefinder_summary }
    
    ch_summary  = ch_summary.mix(serotypefinder_summary)
    ch_versions = ch_versions.mix(SEROTYPEFINDER.out.versions.first())

    log.info "SHIGAPASS is a tool for in silico subtyping of Shigella assemblies, including species assignment, multi-locus sequence typing, and detection of virulence genes."
    SHIGAPASS(ch_ecoli.filter{it})

    SHIGAPASS.out.summary
      .collectFile(name: 'shigapass_summary.csv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigapass")
      .set{ shigapass_summary }


    ch_summary  = ch_summary.mix(shigapass_summary)
    ch_versions = ch_versions.mix(SHIGAPASS.out.versions.first())

  emit:
    for_summary = ch_summary
    versions    = ch_versions

}

if ( ! params.skip_extras ) {
  workflow.onComplete {
    log.info "Subtyping workflow completed at: $workflow.complete"
    log.info "Generated DRPRG summary file: ${params.outdir}/drprg/drprg_summary.tsv"
    log.info "Generated EMMtyper summary file: ${params.outdir}/emmtyper/emmtyper_summary.tsv"
    log.info "Generated KAPTIVE summary file: ${params.outdir}/kaptive/kaptive_summary.txt"
    log.info "Generated Kleborate summary file: ${params.outdir}/kleborate/kleborate_results.tsv"
    log.info "Generated EL GATO summary file: ${params.outdir}/elgato/elgato_summary.tsv"
    log.info "Generated MYKROBE summary file: ${params.outdir}/mykrobe/mykrobe_summary.csv"
    log.info "Generated MENINGOTYPE summary file: ${params.outdir}/meningotype/meningotype_summary.tsv"
    log.info "Generated NGMASTER summary file: ${params.outdir}/ngmaster/ngmaster_summary.tsv"
    log.info "Generated PBPTYPER summary file: ${params.outdir}/pbptyper/pbptyper_summary.tsv"
    log.info "Generated SEQSERO2 summary file: ${params.outdir}/seqsero2/seqsero2_results.txt"
    log.info "Generated SEROTYPEFINDER summary file: ${params.outdir}/serotypefinder/serotypefinder_results.txt"
    log.info "Generated SHIGAPASS summary file: ${params.outdir}/shigapass/*_summary.csv"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
  }
}
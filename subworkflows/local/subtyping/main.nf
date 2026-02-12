
include { DRPRG }          from '../../modules/local/drprg'
include { ELGATO }         from '../../modules/local/elgato'
include { EMMTYPER }       from '../../modules/local/emmtyper'
include { JSON_CONVERT }   from '../../modules/local/local'
include { KAPTIVE }        from '../../modules/local/kaptive'
include { KLEBORATE }      from '../../modules/local/kleborate'
include { MENINGOTYPE }    from '../../modules/local/meningotype'
include { MYKROBE }        from '../../modules/local/mykrobe'
include { NGMASTER }       from '../../modules/local/ngmaster'
include { PBPTYPER }       from '../../modules/local/pbptyper'
include { SEQSERO2 }       from '../../modules/local/seqsero2'
include { SEROTYPEFINDER } from '../../modules/local/serotypefinder'
include { SHIGAPASS }      from '../../modules/local/shigapass'


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


    DRPRG(ch_myco)

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

    EMMTYPER(ch_gas.combine(summfle_script)) 

    EMMTYPER.out.collect
      .collectFile(name: 'emmtyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/emmtyper")
      .set{ emmtyper_summary }

    ch_summary  = ch_summary.mix(emmtyper_summary)
    ch_versions = ch_versions.mix(EMMTYPER.out.versions.first())

    KAPTIVE(ch_vibrio)      

    KAPTIVE.out.collect
      .collectFile(name: 'kaptive_summary.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kaptive")
      .set{ kaptive_summary }
    
    ch_summary  = ch_summary.mix(kaptive_summary)
    ch_versions = ch_versions.mix(KAPTIVE.out.versions.first())

    KLEBORATE(ch_kleb.combine(summfle_script))

    KLEBORATE.out.collect
      .collectFile(name: 'kleborate_results.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kleborate")
      .set{ kleborate_summary }
    
    ch_summary  = ch_summary.mix(kleborate_summary)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())

    ELGATO(ch_legionella)

    ELGATO.out.collect
      .collectFile(name: 'elgato_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/elgato")
      .set{ elgato_summary }

    ch_summary = ch_summary.mix(elgato_summary)
    ch_versions = ch_versions.mix(ELGATO.out.versions.first())

    MYKROBE(ch_myco)

    MYKROBE.out.collect
      .collectFile(name: 'mykrobe_summary.csv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/mykrobe")
      .set{ mykrobe_summary }

    ch_summary  = ch_summary.mix(mykrobe_summary)
    ch_versions = ch_versions.mix(MYKROBE.out.versions.first())

    MENINGOTYPE(ch_gc)

    MENINGOTYPE.out.files
      .collectFile(name: 'meningotype_summary.tsv',
        keepHeader: true,
        sort: {file -> file.text },
        storeDir: "${params.outdir}/meningotype")
      .set{ meningotype_summary }

    ch_summary  = ch_summary.mix(meningotype_summary)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions.first())

    NGMASTER(ch_gc.combine(summfle_script))

    NGMASTER.out.collect
      .collectFile(name: 'ngmaster_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/ngmaster")
      .set{ ngmaster_summary }

    ch_summary  = ch_summary.mix(ngmaster_summary)
    ch_versions = ch_versions.mix(NGMASTER.out.versions.first())

    PBPTYPER(ch_strep)

    PBPTYPER.out.collect
      .collectFile(name: 'pbptyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/pbptyper")
      .set{ pbptyper_summary }

    ch_summary  = ch_summary.mix(pbptyper_summary)
    ch_versions = ch_versions.mix(PBPTYPER.out.versions.first())

    SEQSERO2(ch_salmonella)

    SEQSERO2.out.collect
      .collectFile(name: 'seqsero2_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/seqsero2")
      .set{ seqsero2_summary }

    ch_summary  = ch_summary.mix(seqsero2_summary)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions.first())

    SEROTYPEFINDER(ch_ecoli.combine(summfle_script))

    SEROTYPEFINDER.out.collect
      .collectFile(name: 'serotypefinder_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/serotypefinder")
      .set{ serotypefinder_summary }
    
    ch_summary  = ch_summary.mix(serotypefinder_summary)
    ch_versions = ch_versions.mix(SEROTYPEFINDER.out.versions.first())

    SHIGAPASS(ch_ecoli.combine(summfle_script))

    SHIGAPASS.out.collect
      .collectFile(name: 'shigapass_hits.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigapass")
      .set{ shigapass_hits }

    SHIGAPASS.out.files
      .collectFile(name: 'shigapass_summary.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigapass")
      .set{ shigapass_summary }

    ch_summary  = ch_summary.mix(shigapass_hits).mix(shigapass_summary)
    ch_versions = ch_versions.mix(SHIGAPASS.out.versions.first())

  emit:
    for_summary = ch_summary.collect()
    versions    = ch_versions
}

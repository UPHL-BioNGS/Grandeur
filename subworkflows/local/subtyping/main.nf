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
include { SEQSERO2S }      from '../../../modules/local/seqsero2s'
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
    ch_acinetobacter
    summfle_script
    jsoncon_script

  main:
    ch_summary  = channel.empty()
    ch_versions = channel.empty()

    log.info """

Running subtyping analysis. This workflow will perform in silico subtyping of assemblies 
with a variety of tools, depending on the species of interest.

More information can be found at https://github.com/UPHL-BioNGS/Grandeur/wiki/information.

If a desired sub-typing tool is not included here, please contact the developers or 
submit an issue on GitHub at https://github.com/UPHL-BioNGS/Grandeur/issues

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ DRPRG             ┃ Predicting the drug resistance phenotype of Mycobacterium          ┃
┃                   ┃ tuberculosis from whole genome sequencing data.                    ┃
┃ EMMTYPER          ┃ In silico emm typing of Streptococcus pyogenes assemblies.         ┃
┃ KAPTIVE           ┃ In silico K and O locus typing of Vibrio assemblies, with a focus  ┃
┃                   ┃ on Vibrio parahaemolyticus.                                        ┃
┃ KLEBORATE         ┃ In silico subtyping of Klebsiella assemblies, including species    ┃
┃                   ┃ assignment, multi-locus sequence typing, K and O locus typing, and ┃
┃                   ┃ detection of virulence and AMR genes.                              ┃
┃ EL GATO           ┃ In silico subtyping of Legionella assemblies, including species    ┃
┃                   ┃ assignment, multi-locus sequence typing, and detection of          ┃
┃                   ┃ virulence genes.                                                   ┃
┃ MYKROBE           ┃ Predicting the drug resistance phenotype of Mycobacterium from     ┃
┃                   ┃ whole genome sequencing data.                                      ┃
┃ MENINGOTYPE       ┃ In silico subtyping of Neisseria meningitidis assemblies,          ┃
┃                   ┃ including multi-locus sequence typing and PorA and FetA subtyping. ┃
┃ NGMASTER          ┃ NG-MAST typing of Neisseria gonorrhoeae assemblies.                ┃
┃ PBPTYPER          ┃ In silico penicillin binding protein (PBP) typer of Streptococcus  ┃
┃                   ┃ pneumoniae assemblies.                                             ┃
┃ SEQSERO2          ┃ In silico serotyping of Salmonella assemblies.                     ┃
┃ SEQSERO2S         ┃ In silico serotyping of Salmonella assemblies.                     ┃
┃ SEROTYPEFINDER    ┃ In silico serotyping of Escherichia coli assemblies.               ┃
┃ SHIGAPASS         ┃ In silico subtyping of Shigella assemblies, including species      ┃
┃                   ┃ assignment, multi-locus sequence typing, and detection of          ┃
┃                   ┃ virulence genes.                                                   ┃
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""


    DRPRG(ch_myco.filter{it -> it})

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

    EMMTYPER(ch_gas.filter{it -> it}.combine(summfle_script)) 

    EMMTYPER.out.collect
      .collectFile(name: 'emmtyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/emmtyper")
      .set{ emmtyper_summary }

    ch_summary  = ch_summary.mix(emmtyper_summary)
    ch_versions = ch_versions.mix(EMMTYPER.out.versions.first())

    KAPTIVE(ch_vibrio.mix(ch_acinetobacter).filter{it -> it})      

    KAPTIVE.out.collect
      .collectFile(name: 'kaptive_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kaptive")
      .set{ kaptive_summary }
    
    ch_summary  = ch_summary.mix(kaptive_summary)
    ch_versions = ch_versions.mix(KAPTIVE.out.versions.first())

    KLEBORATE(ch_kleb.mix(ch_ecoli).filter{it -> it}.combine(summfle_script))

    KLEBORATE.out.collect
      .collectFile(name: 'kleborate_results.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kleborate")
      .set{ kleborate_summary }
    
    ch_summary  = ch_summary.mix(kleborate_summary)
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())

    ELGATO(ch_legionella.filter{it -> it})

    ELGATO.out.collect
      .collectFile(name: 'elgato_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/elgato")
      .set{ elgato_summary }

    ch_summary = ch_summary.mix(elgato_summary)
    ch_versions = ch_versions.mix(ELGATO.out.versions.first())

    MYKROBE(ch_myco.filter{it -> it})

    MYKROBE.out.collect
      .collectFile(name: 'mykrobe_summary.csv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/mykrobe")
      .set{ mykrobe_summary }

    ch_summary  = ch_summary.mix(mykrobe_summary)
    ch_versions = ch_versions.mix(MYKROBE.out.versions.first())

    MENINGOTYPE(ch_gc.filter{it -> it})

    MENINGOTYPE.out.summary
      .collectFile(name: 'meningotype_summary.tsv',
        keepHeader: true,
        sort: {file -> file.text },
        storeDir: "${params.outdir}/meningotype")
      .set{ meningotype_summary }

    ch_summary  = ch_summary.mix(meningotype_summary)
    ch_versions = ch_versions.mix(MENINGOTYPE.out.versions.first())

    NGMASTER(ch_gc.filter{it -> it})

    NGMASTER.out.collect
      .collectFile(name: 'ngmaster_summary.csv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/ngmaster")
      .set{ ngmaster_summary }

    ch_summary  = ch_summary.mix(ngmaster_summary)
    ch_versions = ch_versions.mix(NGMASTER.out.versions.first())

    PBPTYPER(ch_strep.filter{it -> it})

    PBPTYPER.out.collect
      .collectFile(name: 'pbptyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/pbptyper")
      .set{ pbptyper_summary }

    ch_summary  = ch_summary.mix(pbptyper_summary)
    ch_versions = ch_versions.mix(PBPTYPER.out.versions.first())

    SEQSERO2(ch_salmonella.filter{it -> it})

    SEQSERO2.out.collect
      .collectFile(name: 'seqsero2_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/seqsero2")
      .set{ seqsero2_summary }

    ch_summary  = ch_summary.mix(seqsero2_summary)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions.first())

    SEQSERO2S(ch_salmonella.filter{it -> it})

    SEQSERO2S.out.collect
      .collectFile(name: 'seqsero2s_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/seqsero2s")
      .set{ seqsero2s_summary }

    ch_summary  = ch_summary.mix(seqsero2s_summary)
    ch_versions = ch_versions.mix(SEQSERO2S.out.versions.first())


    SEROTYPEFINDER(ch_ecoli.filter{it -> it}.combine(summfle_script))

    SEROTYPEFINDER.out.collect
      .collectFile(name: 'serotypefinder_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/serotypefinder")
      .set{ serotypefinder_summary }
    
    ch_summary  = ch_summary.mix(serotypefinder_summary)
    ch_versions = ch_versions.mix(SEROTYPEFINDER.out.versions.first())

    SHIGAPASS(ch_ecoli.filter{it -> it})

    SHIGAPASS.out.summary
      .collectFile(name: 'shigapass_summary.tsv',
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
    log.info """------------------------------------------------------

SUBTYPING subworkflow completed at: $workflow.complete

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Subworkflow Output Files                              ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│   ${params.outdir.padRight(52)}│
│    ├── drprg                                          │
│    │   └── drprg_summary.tsv                          │
│    ├── emmtyper                                       │
│    │   └── emmtyper_summary.tsv                       │
│    ├── kaptive                                        │
│    │   └── kaptive_summary.txt                        │
│    ├── kleborate                                      │
│    │   └── kleborate_results.tsv                      │
│    ├── elgato                                         │
│    │   └── elgato_summary.tsv                         │
│    ├── mykrobe                                        │
│    │   └── mykrobe_summary.csv                        │
│    ├── meningotype                                    │
│    │   └── meningotype_summary.tsv                    │
│    ├── ngmaster                                       │
│    │   └── ngmaster_summary.csv                       │
│    ├── pbptyper                                       │
│    │   └── pbptyper_summary.tsv                       │
│    ├── seqsero2                                       │
│    │   └── seqsero2_results.txt                       │
│    ├── serotypefinder                                 │
│    │   └── serotypefinder_results.txt                 │
│    └── shigapass                                      │
│        └── shigapass_summary.csv                      │
└───────────────────────────────────────────────────────┘

Input files are divided by organism for these processes, 
so not all files will be present.

------------------------------------------------------
"""
  }
}

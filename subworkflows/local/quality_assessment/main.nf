include { AMRFINDER }      from '../../../modules/local/amrfinder'
include { CHECKM2 }         from '../../../modules/local/checkm2'
include { FASTQC }         from '../../../modules/local/fastqc'
include { MLST }           from '../../../modules/local/mlst'
include { PLASMIDFINDER }  from '../../../modules/local/plasmidfinder'
include { QUAST }          from '../../../modules/local/quast'

workflow QUALITY_ASSESSMENT {
    take:
    ch_raw_reads
    ch_clean_reads
    ch_fastas_without_reads
    ch_all_fastas
    ch_reads_contigs
    ch_contigs_org
    ch_checkm2_db
    summfle_script

    main:
    ch_for_multiqc = channel.empty()
    ch_versions    = channel.empty()
    ch_summary     = channel.empty()
    ch_bams        = channel.empty()


    log.info """

Running quality assessment on the reads and assemblies. 

This workflow will perform quality control on the reads with FastQC, but the remaining 
processes of QUAST, CHECKM2, AMRFINDER, and PLASMIDFINDER will be run on generated 
assemblies as well as those specified with an input file designated with 
'params.fasta_list'."

Relevant params and their values:
- 'params.checkm2_db' : ${params.checkm2_db}
    - Set to CHECKM2 database file

┏━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process           ┃ description                                                        ┃
┣━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ FASTQC            ┃ Running quality assessment on the reads with FastQC. This will be  ┃
┃                   ┃ performed on all read files provided, including those specified in ┃
┃                   ┃ the sample sheet, those provided with 'params.reads', and those    ┃
┃                   ┃ downloaded from SRA with 'params.sra_accessions'.                  ┃
┃ AMRFINDER         ┃ AMRFINDERPLUS is a tool for in silico detection of antimicrobial   ┃
┃                   ┃ resistance genes and point mutations in assemblies.                ┃
┃ QUAST             ┃ QUAST is a tool for assessing the quality of genome assemblies by  ┃
┃                   ┃ comparing them to reference genomes and calculating various        ┃
┃                   ┃ assembly metrics.                                                  ┃
┃ MLST              ┃ MLST is a tool for in silico multi-locus sequence typing of        ┃
┃                   ┃ assemblies.                                                        ┃
┃ PLASMIDFINDER     ┃ PLASMIDFINDER is a tool for in silico detection of plasmids in     ┃
┃                   ┃ assemblies.                                                        ┃
┃ CHECKM2           ┃ CHECKM2 is a tool for assessing the quality of genome assemblies   ┃
┃                   ┃ by estimating completeness and contamination based on lineage-     ┃
┃                   ┃ specific marker genes.                                             ┃
┗━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""



    // fastq files only, so hidden if only fasta files are provided
    if ( params.sample_sheet || params.reads || params.sra_accessions ) {
        FASTQC(ch_raw_reads)
        ch_versions    = ch_versions.mix(FASTQC.out.versions.first())
        ch_for_multiqc = ch_for_multiqc.mix(FASTQC.out.for_multiqc)

        FASTQC.out.collect
            .collectFile(name: "fastqc_summary.csv",
                keepHeader: true,
                sort: { file -> file.text },
                storeDir: "${params.outdir}/fastqc")
            .set{ fastqc_summary }

        ch_summary = ch_summary.mix(fastqc_summary)

    }

    AMRFINDER(ch_contigs_org)

    AMRFINDER.out.collect
        .collectFile(name: 'amrfinderplus.txt',
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/amrfinder")
        .set{ amrfinderplus_summary }

    ch_summary  = ch_summary.mix(amrfinderplus_summary)
    ch_versions = ch_versions.mix(AMRFINDER.out.versions.first())

    QUAST(ch_reads_contigs.filter{it})
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    QUAST.out.collect
        .collectFile(name: "quast_report.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/quast")
        .set{ quast_summary }

    ch_summary = ch_summary.mix(quast_summary)
    ch_for_multiqc = ch_for_multiqc.mix(QUAST.out.for_multiqc)

    QUAST.out.collect_contig
        .collectFile(name: "quast_contig_report.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/quast")
        .set{ quast_contig_summary }
    ch_summary = ch_summary.mix(quast_contig_summary)

    

    MLST(ch_all_fastas.combine(summfle_script))
    ch_versions = ch_versions.mix(MLST.out.versions.first())

    MLST.out.collect
        .collectFile(name: "mlst_summary.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/mlst")
        .set{ mlst_summary }
    ch_summary = ch_summary.mix(mlst_summary)

    

    PLASMIDFINDER(ch_all_fastas.combine(summfle_script))
    ch_versions = ch_versions.mix(PLASMIDFINDER.out.versions.first())

    PLASMIDFINDER.out.collect
        .collectFile(name: "plasmidfinder_result.tsv",
            keepHeader: true,
            sort: { file -> file.text },
            storeDir: "${params.outdir}/plasmidfinder")
        .set{ plasmidfinder_summary }
    ch_summary = ch_summary.mix(plasmidfinder_summary)


    if (params.checkm2_db) {    
        CHECKM2(ch_all_fastas.combine(ch_checkm2_db))
        ch_versions = ch_versions.mix(CHECKM2.out.versions.first())

        CHECKM2.out.report
            .collectFile(name: "checkm2_summary.tsv",
                keepHeader: true,
                sort: { file -> file.text },
                storeDir: "${params.outdir}/checkm2")
            .set{ checkm2_summary }
        ch_summary = ch_summary.mix(checkm2_summary)
    } 

    emit:
    for_summary = ch_summary
    for_multiqc = ch_for_multiqc
    versions    = ch_versions
}

if ( ! params.skip_extras ) {
    workflow.onComplete {
        log.info """------------------------------------------------------

QUALITY ASSESSMENT subworkflow completed at: $workflow.complete

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Subworkflow Output Files                              ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│   'params.outdir'                                     │
│    ├── fastqc                                         │
│    │   └── fastqc_summary.csv                         │"""
        if ( params.checkm2_db ) {
            log.info """│    ├── checkm2                                        │
│    │   └── checkm2_summary.tsv                        │"""
        }

        log.info """│    ├── mlst                                           │
│    │   └── mlst_summary.tsv                           │
│    ├── plasmidfinder                                  │
│    │   └── plasmidfinder_result.tsv                   │
│    └── quast                                          │
│        └── quast_report.tsv                           │
└───────────────────────────────────────────────────────┘

------------------------------------------------------
"""
    }
}

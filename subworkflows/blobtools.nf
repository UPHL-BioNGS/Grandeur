include { bwa }                     from '../modules/bwa'        addParams(bwa_options: params.bwa_options )
include { sort }                    from '../modules/samtools'   addParams(samtools_sort_options: params.samtools_sort_options )
include { blastn }                  from '../modules/blast'      addParams(local_db_type: params.local_db_type )
include { blobtools_create; blobtools_view; blobtools_blobtools } from '../modules/blobtools'  addParams(blobtools_create_options: params.blobtools_create_options, blobtools_view_options: params.blobtools_view_options, blobtools_plot_options: params.blobtools_plot_options)

workflow blobtools {
  take:
    clean_reads
    contigs
    blast_db
  main:
    bwa(clean_reads.join(contigs, by: 0))
    blastn(contigs.combine(blast_db))
    sort(bwa.out.sam)
    blobtools_create(contigs.join(blastn.out.blastn, by: 0).join(sort.out.bam, by: 0))
    blobtools_view(blobtools_create.out.json)
    blobtools_blobtools(blobtools_create.out.json)
  emit:
    species = blobtools_blobtools.out.species
    perc    = blobtools_blobtools.out.perc
}

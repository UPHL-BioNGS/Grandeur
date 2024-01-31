include { download_sra }  from '../modules/local/download_sra' addParams(params)

workflow test {
    take:
        ch_accessions
  
    main:
        download_sra(ch_accessions)
    
    emit:
        fastq    = download_sra.out.fastq
        versions = download_sra.out.versions.first()
}
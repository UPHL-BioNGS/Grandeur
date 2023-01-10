include { download_sra }  from '../modules/download_sra' addParams(params)

workflow test {
    take:
        ch_accessions
  
    main:
        download_sra(ch_accessions)
    
    emit:
        fastq = download_sra.out.fastq

}
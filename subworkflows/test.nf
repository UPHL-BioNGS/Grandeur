include { download_sra }  from '../modules/local/local' addParams(params)

workflow test {
    take:
        ch_accessions
  
    main:
        download_sra(ch_accessions)

        download_sra.out.fastq
            .map { it ->
                meta = [id:it[0]] 
                tuple( meta, [file(it[1][0]), file(it[1][1])])}
            .set { ch_fastq }
    
    emit:
        fastq    = ch_fastq
}
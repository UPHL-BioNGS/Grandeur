include { seqyclean }                                from '../modules/seqyclean'
include { spades }                               from '../modules/spades'

workflow de_novo_alignment {
  take: reads
  main:
    seqyclean(reads)
    spades(clean_reads)
}

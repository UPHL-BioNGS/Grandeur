
import groovy.json.JsonSlurper

include { AMRFINDER }      from '../../modules/local/amrfinderplus'
include { DRPRG }          from '../../modules/local/drprg'
include { ELGATO }         from '../../modules/local/elgato'
include { EMMTYPER }       from '../../modules/local/emmtyper'
include { KAPTIVE }        from '../../modules/local/kaptive'
include { KLEBORATE }      from '../../modules/local/kleborate'
include { MENINGOTYPE }    from '../../modules/local/meningotype'
include { MYKROBE }        from '../../modules/local/mykrobe'
include { PBPTYPER }       from '../../modules/local/pbptyper'
include { SEQSERO2 }       from '../../modules/local/seqsero2'
include { SEROTYPEFINDER } from '../../modules/local/serotypefinder'
include { SHIGATYPER }     from '../../modules/local/shigatyper'

def flagOrg(org_files, phrases) {
  def found = false
  org_files.each { org_file ->
    if (org_file && org_file.exists()) {
      def count = 0
      org_file.withReader { reader ->
        while(reader.ready() && count < 10 && !found) {
          def line = reader.readLine()
          count = count + 1
          phrases.each { phrase ->
            if (line.toString().contains(phrase)) {
              if (org_file.getName().contains('fastani')) {
                def columns = line.split(',')
                if (columns.size() >= 4 && columns[3].toFloat() > 90) {
                  found = true
                }
              } else {
                found = true
              }
            }     
          }
        }
      }
    }
  }
  return found
}

def topOrg(org_files) {
  def fastani_file = org_files[0]
  def blobtools_file = org_files[1]
  def kraken2_file = org_files[2]
  def mash_file = org_files[3]
  
  def genus   = 'unknown'
  def species = 'unknown'

  if (fastani_file && fastani_file.exists() && genus == 'unknown') {
    def lines = fastani_file.readLines()
    if (lines.size() > 1) {
      def secondLine = lines[1].split(',')
      if (secondLine.size() >= 4 && secondLine[3].toFloat() > 90) {
        def hit = secondLine[2].trim()+ '_unknown'
        genus   = hit.split('_')[0]
        species = hit.split('_')[1]
      }
    }
  } 
  
  if (blobtools_file && blobtools_file.exists() && genus == 'unknown') {
    blobtools_file.withReader { reader -> 
      while(reader.ready()) {
        def line = reader.readLine()
        def columns = line.split('\t')
        if (columns.size() > 1) {
          def hit = columns[1].trim()
          if (!['name', 'all', 'no-hit', 'undel'].contains(hit)) {
            if (columns.size() == 14 && columns[-1].toFloat() > 50) {
              def name = hit + '_unknown'
              genus   = name.split('_')[0]
              species = name.split('_')[1]
              return
            }
          }
        }
      }
    }
  }
  
  if (kraken2_file && kraken2_file.exists() && genus == 'unknown') {
    def lines = kraken2_file.readLines()
    if (lines.size() > 1) {
      def secondLine = lines[1].split(',')
      if (secondLine.size() >= 2 && secondLine[1].toFloat() > 50) {
        def hit = secondLine[-1].trim() + '_unknown'
        genus   = hit.split('_')[0]
        species = hit.split('_')[1]
      }
    }
  } 
  
  if (mash_file && mash_file.exists() && genus == 'unknown') {
    def lines = mash_file.readLines()
    if (lines.size() > 1) {
      def secondLine = lines[1].split(',')
      if (secondLine.size() >= 4 && secondLine[3].toFloat() < 0.1) {
        def hit = secondLine[-1].trim() + '_unknown'
        genus   = hit.split('_')[0]
        species = hit.split('_')[1]
      }
    }
  }
  return [genus, species]
}


// def jsonConvert(it) {


//   def filename = it[1]

//   // Read and parse the JSON file
//   def jsonSlurper = new JsonSlurper()
//   def contents = [:]

//   new File(filename).withReader { reader ->
//       contents = jsonSlurper.parse(reader)
//   }

//   // Extract data from JSON
//   def sample = contents.sample
//   def present = contents.genes.present.join(",")
//   def absent = contents.genes.absent.join(",")
//   def susceptibility = contents.susceptibility.keySet().join(",")

//   // Create the summary string
//   def summary = "sample\tgenes_present\tgenes_absent\tsusceptibility\n"
//   summary += "${sample}\t${present}\t${absent}\t${susceptibility}"

//   // Return the summary string
//   println summary
// }

workflow INFO {
  take:
    ch_contigs
    ch_flag
    summfle_script
    jsoncon_script

  main:
    ch_summary  = Channel.empty()
    ch_versions = Channel.empty()

    // species specific
    // branch + join = faster than groupTuple
    ch_flag
      .branch {
        blobtools: it[1] =~ /blobtools.txt/
        kraken2:   it[1] =~ /kraken2.csv/
        mash:      it[1] =~ /mash.csv/
        fastani:   it[1] =~ /fastani.csv/
      }
      .set { ch_flag_branch }
  

    ch_contigs
      .filter{it[1] != null}
      .join(ch_flag_branch.fastani,   by:0, failOnMismatch: false, remainder: true)
      .join(ch_flag_branch.blobtools, by:0, failOnMismatch: false, remainder: true)
      .join(ch_flag_branch.kraken2,   by:0, failOnMismatch: false, remainder: true)
      .join(ch_flag_branch.mash,      by:0, failOnMismatch: false, remainder: true)
      .filter{it[1] != null}
      .map{ it -> tuple(it[0], it[1], [it[2], it[3], it[4], it[5]])}
      .set {ch_for_flag}

    // Yersinia
    ch_for_flag
      .filter{flagOrg(it[2], ['Yersinia'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_yersinia}

    // Salmonella
    ch_for_flag
      .filter{flagOrg(it[2], ['Salmonella'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_salmonella}

    // E. coli and Shigella
    ch_for_flag
      .filter{flagOrg(it[2], ['Escherichia', 'Shigella'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_ecoli}

    // Klebsiella
    ch_for_flag
      .filter{flagOrg(it[2], ['Klebsiella', 'Enterobacter', 'Serratia'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_kleb}

    // Group A Strep
    ch_for_flag
      .filter{flagOrg(it[2], ['Streptococcus'])}
      .filter{flagOrg(it[2], ['pyogenes', 'dysgalactiae', 'anginosus'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_gas}

    // Streptococcus pneumoniae
    ch_for_flag
      .filter{flagOrg(it[2], ['Streptococcus'])}
      .filter{flagOrg(it[2], ['pneumoniae'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_strep}

    // Legionella
    ch_for_flag
      .filter{flagOrg(it[2], ['Legionella'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_legionella}

    // Vibrio
    ch_for_flag
      .filter{flagOrg(it[2], ['Vibrio'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_vibrio}

    // Acinetobacter
    ch_for_flag
      .filter{flagOrg(it[2], ['Acinetobacter'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_acinetobacter}

    // Mycobacterium/Mycobacteria
    ch_for_flag
      .filter{flagOrg(it[2], ['Mycobacteri'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_myco}


    // Neisseria meningitidis
    ch_for_flag
      .filter{flagOrg(it[2], ['Neisseria'])}
      .map { it -> tuple(it[0], it[1])}
      .set {ch_gc}

    // Getting the top organism for each sample
    // for amrfinderplus
    // for prokka
    ch_for_flag
      .map {
        it ->
          def genus_species = topOrg(it[2])
          tuple (it[0], it[1], genus_species[0], genus_species[1])
      }
      .set { ch_organism }

    AMRFINDER(ch_organism)

    AMRFINDER.out.collect
      .collectFile(name: 'amrfinderplus.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/ncbi-AMRFinderplus")
      .set{ amrfinderplus_summary }

    ch_summary  = ch_summary.mix(amrfinderplus_summary)
    ch_versions = ch_versions.mix(AMRFINDER.out.versions.first())

    DRPRG(ch_myco)

    // DRPRG.out.json
    //   .map { it -> jsonConvert(it) }
    //   .collectFile(name: 'drprg_summary.tsv',
    //     keepHeader: true,
    //     sort: { file -> file.text },
    //     storeDir: "${params.outdir}/drprg")
    //   .set{ drprg_summary }

    //ch_summary  = ch_summary.mix(drprg_summary)
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

    SHIGATYPER(ch_ecoli.combine(summfle_script))

    SHIGATYPER.out.collect
      .collectFile(name: 'shigatyper_hits.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigatyper")
      .set{ shigatyper_hits }

    SHIGATYPER.out.files
      .collectFile(name: 'shigatyper_summary.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigatyper")
      .set{ shigatyper_summary }

    ch_summary  = ch_summary.mix(shigatyper_hits).mix(shigatyper_summary)
    ch_versions = ch_versions.mix(SHIGATYPER.out.versions.first())

  emit:
    for_summary = ch_summary.collect()
    versions    = ch_versions
}

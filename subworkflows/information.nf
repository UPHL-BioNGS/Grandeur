include { amrfinderplus }  from '../modules/local/amrfinderplus'  addParams(params)
include { drprg }          from '../modules/local/drprg'          addParams(params)
include { elgato }         from '../modules/local/elgato'         addParams(params)
include { emmtyper }       from '../modules/local/emmtyper'       addParams(params)
include { json_convert }   from '../modules/local/local'          addParams(params)
include { kaptive }        from '../modules/local/kaptive'        addParams(params)
include { kleborate }      from '../modules/local/kleborate'      addParams(params)
include { mykrobe }        from '../modules/local/mykrobe'        addParams(params)
include { pbptyper }       from '../modules/local/pbptyper'       addParams(params)
include { seqsero2 }       from '../modules/local/seqsero2'       addParams(params)
include { serotypefinder } from '../modules/local/serotypefinder' addParams(params)
include { shigatyper }     from '../modules/local/shigatyper'     addParams(params)

def flagOrg(org_files, phrases) {
  def found = false
  org_files.each { org_file ->
    if (org_file && org_file.exists()) {
      def count = 0
      org_file.withReader { reader ->
        while (reader.ready() && count < 10 && !found) {
          def line = reader.readLine()
          count++
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
        hit     = secondLine[2].trim()+ '_unknown'
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
          hit = columns[1].trim()
          if (!['name', 'all', 'no-hit', 'undel'].contains(hit)) {
            if (columns.size() == 14 && columns[-1].toFloat() > 50) {
              name    = hit + '_unknown'
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
        hit     = secondLine[-1].trim() + '_unknown'
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
        hit     = secondLine[-1].trim() + '_unknown'
        genus   = hit.split('_')[0]
        species = hit.split('_')[1]
      }
    }
  }
  return [genus, species]
}

workflow information {
  take:
    ch_contigs
    ch_flag
    summfle_script
    jsoncon_script

  main:

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

    // Getting the top organism for each sample
    // for amrfinderplus
    // for prokka
    ch_for_flag
      .map {
        it ->
          genus_species = topOrg(it[2])
          tuple (it[0], it[1], genus_species[0], genus_species[1])
      }
      .set { ch_organism }

    amrfinderplus(ch_organism)
    drprg(ch_myco)
    emmtyper(ch_gas.combine(summfle_script)) 
    kaptive(ch_vibrio)      
    kleborate(ch_kleb.combine(summfle_script))
    elgato(ch_legionella)
    mykrobe(ch_myco)
    pbptyper(ch_strep)
    seqsero2(ch_salmonella)
    serotypefinder(ch_ecoli.combine(summfle_script))
    shigatyper(ch_ecoli.combine(summfle_script))

    json_convert(drprg.out.json.combine(jsoncon_script))

    amrfinderplus.out.collect
      .collectFile(name: 'amrfinderplus.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/ncbi-AMRFinderplus")
      .set{ amrfinderplus_summary }

    json_convert.out.collect
      .filter( ~/.*drprg.tsv/ )
      .collectFile(name: 'drprg_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/drprg")
      .set{ drprg_summary }

    elgato.out.collect
      .collectFile(name: 'elgato_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/elgato")
      .set{ elgato_summary }

    emmtyper.out.collect
      .collectFile(name: 'emmtyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/emmtyper")
      .set{ emmtyper_summary }

    kaptive.out.collect
      .collectFile(name: 'kaptive_summary.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kaptive")
      .set{ kaptive_summary }
    
    kleborate.out.collect
      .collectFile(name: 'kleborate_results.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/kleborate")
      .set{ kleborate_summary }

    mykrobe.out.collect
      .collectFile(name: 'mykrobe_summary.csv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/mykrobe")
      .set{ mykrobe_summary }

    pbptyper.out.collect
      .collectFile(name: 'pbptyper_summary.tsv',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/pbptyper")
      .set{ pbptyper_summary }

    seqsero2.out.collect
      .collectFile(name: 'seqsero2_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/seqsero2")
      .set{ seqsero2_summary }

    serotypefinder.out.collect
      .collectFile(name: 'serotypefinder_results.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/serotypefinder")
      .set{ serotypefinder_summary }

    shigatyper.out.collect
      .collectFile(name: 'shigatyper_hits.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigatyper")
      .set{ shigatyper_hits }

    shigatyper.out.files
      .collectFile(name: 'shigatyper_summary.txt',
        keepHeader: true,
        sort: { file -> file.text },
        storeDir: "${params.outdir}/shigatyper")
      .set{ shigatyper_summary }

    amrfinderplus_summary
      .mix(drprg_summary)
      .mix(elgato_summary)
      .mix(emmtyper_summary)
      .mix(kaptive_summary)
      .mix(kleborate_summary)
      .mix(mykrobe_summary)
      .mix(pbptyper_summary)
      .mix(seqsero2_summary)
      .mix(serotypefinder_summary)
      .mix(shigatyper_hits)
      .mix(shigatyper_summary)
      .set { for_summary }

    amrfinderplus.out.versions.first()
      .mix(drprg.out.versions)
      .mix(elgato.out.versions)
      .mix(emmtyper.out.versions)
      .mix(kaptive.out.versions)
      .mix(kleborate.out.versions)
      .mix(mykrobe.out.versions)
      .mix(pbptyper.out.versions)
      .mix(seqsero2.out.versions)
      .mix(serotypefinder.out.versions)
      .mix(shigatyper.out.versions)
      .set { for_versions }

  emit:
    for_summary = for_summary.collect()
    versions    = for_versions
}

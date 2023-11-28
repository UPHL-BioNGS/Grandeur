#!/usr/bin/env nextflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Welcome to this workflow! Issues and contributions are gladly accepted at https://github.com/UPHL-BioNGS/Grandeur .

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

println("Currently using the Grandeur workflow for use with microbial sequencing. The view is great from 8299 feet (2530 meters) above sea level.\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: ${workflow.manifest.version}")
println("")

nextflow.enable.dsl               = 2

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Getting config file

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.config_file                = false
if ( params.config_file ) {
  def src = new File("${workflow.projectDir}/configs/grandeur_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}

params.fastqscan_options = false
if ( params.fastqscan_options ) {
  println("WARNING : params.fastqscan_options is depreciated" )
}

params.genome_references = false
if (params.genome_references ) {
  println("WARNING : params.genome_references is depreciated" )
}

params.roary_min_genes = false
if (params.genome_references ) {
  println("WARNING : params.roary_min_genes is depreciated, use params.min_core_genes instead!")
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Defining params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

params.outdir                     = "grandeur"
params.maxcpus                    = 12
params.medcpus                    = 4
params.minimum_reads              = 10000

// connecting to Donut Falls (in development)
params.donut_falls_wf             = false

// input files
params.reads                      = workflow.launchDir + "/reads"
params.fastas                     = workflow.launchDir + "/fastas"
params.gff                        = workflow.launchDir + "/gff"
params.sample_sheet               = ""
params.fasta_list                 = ""

// external files
params.kraken2_db                 = ""
params.blast_db                   = ""
params.mash_db                    = ""
params.fastani_ref                = ""
params.fastani_ref_list           = ""
params.genome_sizes               = workflow.projectDir + "/assets/genome_sizes.json"

// for testing
params.sra_accessions             = []

// tool-specific command line options
params.amrfinderplus_options      = ""
params.bbduk_options              = "k=31 hdist=1"
params.bbmap_options              = ""
params.blast_db_type              = "nt"
params.blastn_options             = "-max_target_seqs 10 -max_hsps 1 -evalue 1e-25"
params.blobtools_create_options   = ""
params.blobtools_view_options     = ""
params.blobtools_plot_options     = "--format png -r species"
params.blobtools_bbmap_options    = ""
params.current_datasets           = false
params.datasets_max_genomes       = 5
params.drprg_options              = ""
params.emmtyper_options           = ''
params.extras                     = true
params.fastani_include            = true
params.fastani_options            = ""
params.fasterqdump_options        = ""
params.fastp_options              = "--detect_adapter_for_pe"
params.fastqc_options             = ""
params.heatcluster_options        = "-t png"
params.iqtree2_options            = "-t RANDOM -m GTR+F+I -bb 1000 -alrt 1000"
params.iqtree2_outgroup           = ""
params.kaptive_options            = ''
params.kleborate_options          = "-all"
params.kraken2_options            = ""
params.legsta_options             = ""
params.mash_sketch_options        = "-m 2"
params.mash_dist_options          = "-v 0 -d 0.5"
params.mash_max_hits              = 25
params.mashtree_options           = ""
params.min_core_genes             = 1500
params.mlst_options               = ""
params.msa                        = false
params.multiqc_options            = ""
params.mykrobe_options            = ""
params.panaroo_options            = "--clean-mode strict --remove-invalid-genes"
params.pbptyper_options           = ''
params.plasmidfinder_options      = ""
params.phytreeviz_options         = ""
params.prokka_options             = "--mincontiglen 500 --compliant --locustag locus_tag --centre STAPHB"
params.quast_options              = ""
params.roary_options              = ""
params.seqsero2_options           = "-m a -b mem"
params.serotypefinder_options     = ""
params.shigatyper_options         = ""
params.snp_dists_options          = "-c"
params.spades_options             = "--isolate"

// if ( params.phoenix_wf ) {
//   println "cp ${workflow.projectDir}/../phoenix/bin/* ${workflow.projectDir}/bin/."
//   command = ["sh", "-c", "cp ${workflow.projectDir}/../phoenix/bin/* ${workflow.projectDir}/bin/."]
//   Runtime.getRuntime().exec((String[]) command.toArray())
//   command = ["sh", "-c", "cp ${workflow.projectDir}/../phoenix/lib/* ${workflow.projectDir}/lib/."]
//   Runtime.getRuntime().exec((String[]) command.toArray())
//   include { PHOENIX_EXTERNAL } from workflow.projectDir + "/../phoenix/main" addParams(input: params.input)
//   For my future self, this doesn't work because of some of the nf-core scripts
// }

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Sharing params with subworkflows

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

include { average_nucleotide_identity }   from "./subworkflows/average_nucleotide_identity"   addParams(params)
include { blobtools }                     from "./subworkflows/blobtools"                     addParams(params)
include { de_novo_alignment }             from "./subworkflows/de_novo_alignment"             addParams(params)
include { information }                   from "./subworkflows/information"                   addParams(params)
include { kmer_taxonomic_classification } from "./subworkflows/kmer_taxonomic_classification" addParams(params)
include { min_hash_distance }             from "./subworkflows/min_hash_distance"             addParams(params)
include { phylogenetic_analysis }         from "./subworkflows/phylogenetic_analysis"         addParams(params)
include { report }                        from "./subworkflows/report"                        addParams(params)
include { test }                          from "./subworkflows/test"                          addParams(params)

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for scripts

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Creating the summary files
dataset_script = Channel.fromPath(workflow.projectDir + "/bin/datasets_download.py", type: "file")
evaluat_script = Channel.fromPath(workflow.projectDir + "/bin/evaluate.py",          type: "file")
summary_script = Channel.fromPath(workflow.projectDir + "/bin/summary.py",           type: "file")
summfle_script = Channel.fromPath(workflow.projectDir + "/bin/summary_file.py",      type: "file")

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Included genomes for fastani

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

included_genomes = [
  'Acinetobacter_baumannii_GCF_008632635.1.fna.gz',
  'Acinetobacter_bereziniae_GCF_016576965.1.fna.gz',
  'Acinetobacter_calcoaceticus_GCF_002055515.1.fna.gz',
  'Acinetobacter_guillouiae_GCF_002370525.2.fna.gz',
  'Acinetobacter_haemolyticus_GCF_003323815.1.fna.gz',
  'Acinetobacter_junii_GCF_018336855.1.fna.gz',
  'Acinetobacter_lactucae_GCF_013122135.1.fna.gz',
  'Acinetobacter_nosocomialis_M2_GCF_005281455.1.fna.gz',
  'Acinetobacter_pittii_PHEA-2_GCF_000191145.1.fna.gz',
  'Acinetobacter_proteolyticus_GCF_000367945.1.fna.gz',
  'Acinetobacter_radioresistens_GCF_003258335.1.fna.gz',
  'Acinetobacter_schindleri_GCF_010918895.1.fna.gz',
  'Acinetobacter_seifertii_GCF_016064815.1.fna.gz',
  'Acinetobacter_variabilis_GCF_018409485.1.fna.gz',
  'Alcaligenes_faecalis_GCF_000967305.2.fna.gz',
  'Campylobacter_coli_GCA_008011635.1.fna.gz',
  'Campylobacter_coli_GCF_009730395.1.fna.gz',
  'Campylobacter_fetus_GCF_011600945.2.fna.gz',
  'Campylobacter_fetus_subsp._fetus_82-40_GCF_000015085.1.fna.gz',
  'Campylobacter_fetus_subsp._testudinum_03-427_GCF_000495505.1.fna.gz',
  'Campylobacter_fetus_subsp._venerealis_GCF_000759515.1.fna.gz',
  'Campylobacter_hyointestinalis_subsp._hyointestinalis_LMG_9260_GCF_001643955.1.fna.gz',
  'Campylobacter_hyointestinalis_subsp._lawsonii_GCF_013372165.1.fna.gz',
  'Campylobacter_jejuni_subsp._doylei_269.97_GCF_000017485.1.fna.gz',
  'Campylobacter_jejuni_subsp._jejuni_GCA_008011525.1.fna.gz',
  'Campylobacter_jejuni_subsp._jejuni_NCTC_11168_ATCC_700819_GCF_000009085.1.fna.gz',
  'Campylobacter_lari_RM2100_GCF_000019205.1.fna.gz',
  'Campylobacter_lari_subsp._concheus_LMG_11760_GCF_000816225.1.fna.gz',
  'Campylobacter_peloridis_GCF_014931075.1.fna.gz',
  'Campylobacter_subantarcticus_LMG_24377_GCF_000816305.1.fna.gz',
  'Campylobacter_upsaliensis_GCA_008011615.1.fna.gz',
  'Campylobacter_upsaliensis_GCF_916098265.1.fna.gz',
  'Citrobacter_amalonaticus_GCF_001558935.2.fna.gz',
  'Citrobacter_braakii_GCF_009648935.1.fna.gz',
  'Citrobacter_farmeri_GCF_003938205.1.fna.gz',
  'Citrobacter_freundii_GCF_003812345.1.fna.gz',
  'Citrobacter_gillenii_GCF_003429605.1.fna.gz',
  'Citrobacter_koseri_ATCC_BAA-895_GCF_000018045.1.fna.gz',
  'Citrobacter_murliniae_GCF_004801125.1.fna.gz',
  'Citrobacter_pasteurii_GCF_019047765.1.fna.gz',
  'Citrobacter_portucalensis_GCF_008693605.1.fna.gz',
  'Citrobacter_sedlakii_GCF_018128425.1.fna.gz',
  'Citrobacter_werkmanii_GCF_008693645.1.fna.gz',
  'Citrobacter_youngae_GCF_900638065.1.fna.gz',
  'Clostridium_botulinum_A_str._ATCC_3502_GCF_000063585.1.fna.gz',
  'Clostridium_perfringens_GCF_020138775.1.fna.gz',
  'Clostridium_sporogenes_GCF_020450145.1.fna.gz',
  'Cronobacter_dublinensis_subsp._dublinensis_LMG_23823_GCF_001277235.1.fna.gz',
  'Cronobacter_malonaticus_LMG_23826_GCF_001277215.2.fna.gz',
  'Cronobacter_muytjensii_ATCC_51329_GCF_001277195.1.fna.gz',
  'Cronobacter_sakazakii_GCF_003516125.1.fna.gz',
  'Cronobacter_turicensis_GCF_011605535.1.fna.gz',
  'Cronobacter_universalis_NCTC_9529_GCF_001277175.1.fna.gz',
  'Elizabethkingia_anophelis_R26_GCF_002023665.2.fna.gz',
  'Elizabethkingia_meningoseptica_GCF_002022145.1.fna.gz',
  'Enterobacter_asburiae_GCF_007035805.1.fna.gz',
  'Enterobacter_bugandensis_GCF_015137655.1.fna.gz',
  'Enterobacter_bugandensis_GCF_020042625.1.fna.gz',
  'Enterobacter_cancerogenus_GCF_019665745.1.fna.gz',
  'Enterobacter_chengduensis_GCF_001984825.2.fna.gz',
  'Enterobacter_chuandaensis_GCF_003594915.1.fna.gz',
  'Enterobacter_cloacae_GCF_023702375.1.fna.gz',
  'Enterobacter_cloacae_GCF_905331265.2.fna.gz',
  'Enterobacter_dykesii_GCF_018597265.1.fna.gz',
  'Enterobacter_hormaechei_GCF_019048625.1.fna.gz',
  'Enterobacter_huaxiensis_GCF_003594935.2.fna.gz',
  'Enterobacter_kobei_GCF_000534275.1.fna.gz',
  'Enterobacter_ludwigii_GCF_001750725.1.fna.gz',
  'Enterobacter_mori_GCF_022014715.1.fna.gz',
  'Enterobacter_oligotrophicus_GCF_009176645.1.fna.gz',
  'Enterobacter_quasihormaechei_GCF_004331385.1.fna.gz',
  'Enterobacter_quasimori_GCF_018597345.1.fna.gz',
  'Enterobacter_quasiroggenkampii_GCF_003964805.1.fna.gz',
  'Enterobacter_roggenkampii_GCF_001729805.1.fna.gz',
  'Enterobacter_sichuanensis_GCF_009036245.1.fna.gz',
  'Enterobacter_soli_GCF_000224675.1.fna.gz',
  'Enterobacter_vonholyi_GCF_008364555.1.fna.gz',
  'Enterobacter_wuhouensis_GCF_004331265.1.fna.gz',
  'Enterococcus_faecalis_EnGen0336_GCF_000393015.1.fna.gz',
  'Enterococcus_faecium_GCF_009734005.1.fna.gz',
  'Escherichia_albertii_GCF_016904755.1.fna.gz',
  'Escherichia_albertii_KF1_GCF_000512125.1.fna.gz',
  'Escherichia_coli_O157:H7_str._Sakai_GCF_000008865.2.fna.gz',
  'Escherichia_coli_O27:H7_GCF_002741475.1.fna.gz',
  'Escherichia_coli_str._K-12_substr._MG1655_GCF_000005845.2.fna.gz',
  'Escherichia_fergusonii_ATCC_35469_GCF_000026225.1.fna.gz',
  'Escherichia_fergusonii_GCF_020097475.1.fna.gz',
  'Grimontia_hollisae_GCF_009665295.1.fna.gz',
  'Haemophilus_aegyptius_GCF_900475885.1.fna.gz',
  'Haemophilus_haemolyticus_GCF_900477945.1.fna.gz',
  'Haemophilus_influenzae_GCF_000931575.1.fna.gz',
  'Haemophilus_parainfluenzae_ATCC_33392_GCF_000191405.1.fna.gz',
  'Hafnia_alvei_GCF_011617105.1.fna.gz',
  'Hafnia_paralvei_GCF_020150375.1.fna.gz',
  'Klebsiella_aerogenes_GCF_007632255.1.fna.gz',
  'Klebsiella_michiganensis_GCF_015139575.1.fna.gz',
  'Klebsiella_oxytoca_GCF_003812925.1.fna.gz',
  'Klebsiella_pasteurii_GCF_018139045.1.fna.gz',
  'Klebsiella_pneumoniae_subsp._pneumoniae_HS11286_GCF_000240185.1.fna.gz',
  'Klebsiella_quasipneumoniae_GCF_016415705.1.fna.gz',
  'Klebsiella_variicola_GCF_009648975.1.fna.gz',
  'Kluyvera_ascorbata_GCF_023195735.1.fna.gz',
  'Leclercia_adecarboxylata_GCF_001518835.1.fna.gz',
  'Legionella_pneumophila_GCF_001753085.1.fna.gz',
  'Lelliottia_amnigena_GCF_019355955.1.fna.gz',
  'Listeria_innocua_GCF_009648575.1.fna.gz',
  'Listeria_innocua_GCF_017363615.1.fna.gz',
  'Listeria_innocua_GCF_017363655.1.fna.gz',
  'Listeria_ivanovii_subsp._ivanovii_PAM_55_GCF_000252975.1.fna.gz',
  'Listeria_marthii_GCF_017363645.1.fna.gz',
  'Listeria_monocytogenes_EGD-e_GCF_000196035.1.fna.gz',
  'Listeria_monocytogenes_GCF_001466295.1.fna.gz',
  'Listeria_monocytogenes_GCF_013625895.1.fna.gz',
  'Listeria_monocytogenes_GCF_013625995.1.fna.gz',
  'Listeria_monocytogenes_GCF_013626145.1.fna.gz',
  'Listeria_monocytogenes_GCF_014526935.1.fna.gz',
  'Listeria_seeligeri_GCF_017363605.1.fna.gz',
  'Listeria_welshimeri_GCF_002489005.1.fna.gz',
  'Mixta_calida_GCA_007681265.1.fna.gz',
  'Morganella_morganii_GCF_902387845.1.fna.gz',
  'Mycobacterium_avium_subsp_hominissuis_GCF_022175585.1.fna.gz',
  'Mycobacterium_leprae_GCF_003253775.1.fna.gz',
  'Mycobacterium_marinum_E11_GCF_000723425.2.fna.gz',
  'Mycobacterium_tuberculosis_H37Rv_GCF_000195955.2.fna.gz',
  'Mycobacterium_ulcerans_GCF_020616615.1.fna.gz',
  'Neisseria_gonorrhoeae_GCF_013030075.1.fna.gz',
  'Pantoea_ananatis_PA13_GCF_000233595.1.fna.gz',
  'Photobacterium_damselae_GCF_009665375.1.fna.gz',
  'Pluralibacter_gergoviae_GCF_003019925.1.fna.gz',
  'Proteus_hauseri_GCF_004116975.1.fna.gz',
  'Proteus_mirabilis_HI4320_GCF_000069965.1.fna.gz',
  'Proteus_vulgaris_GCF_000754995.1.fna.gz',
  'Providencia_huaxiensis_GCF_002843235.3.fna.gz',
  'Providencia_rettgeri_GCF_003204135.1.fna.gz',
  'Providencia_stuartii_GCF_023547145.1.fna.gz',
  'Providencia_stuartii_GCF_029277985.1.fna.gz',
  'Pseudescherichia_vulneris_GCF_902164725.1.fna.gz',
  'Pseudomonas_aeruginosa_PAO1_GCF_000006765.1.fna.gz',
  'Pseudomonas_alcaligenes_GCF_001597285.1.fna.gz',
  'Pseudomonas_fluorescens_GCF_900215245.1.fna.gz',
  'Pseudomonas_fulva_GCF_001186195.1.fna.gz',
  'Pseudomonas_furukawaii_GCF_002355475.1.fna.gz',
  'Pseudomonas_multiresinivorans_GCF_012971725.1.fna.gz',
  'Pseudomonas_nitroreducens_GCF_012986205.1.fna.gz',
  'Pseudomonas_otitidis_GCF_011397855.1.fna.gz',
  'Pseudomonas_putida_NBRC_14164_GCF_000412675.1.fna.gz',
  'Ralstonia_pickettii_GCF_902374465.1.fna.gz',
  'Raoultella_ornithinolytica_GCF_901421005.1.fna.gz',
  'Raoultella_planticola_GCF_022637595.1.fna.gz',
  'Salmonella_bongori_N268-08_GCF_000439255.1.fna.gz',
  'Salmonella_enterica_GCA_011388235.1.fna.gz',
  'Salmonella_enterica_subsp._enterica_serovar_Typhimurium_str._LT2_GCF_000006945.2.fna.gz',
  'Salmonella_enterica_subsp._houtenae_GCA_013588055.1.fna.gz',
  'Serratia_marcescens_GCF_003516165.1.fna.gz',
  'Serratia_nematodiphila_GCF_004768745.1.fna.gz',
  'Shigella_boydii_GCF_002290485.1.fna.gz',
  'Shigella_dysenteriae_GCF_022354085.1.fna.gz',
  'Shigella_flexneri_2a_str._301_GCF_000006925.2.fna.gz',
  'Shigella_sonnei_GCF_013374815.1.fna.gz',
  'Staphylococcus_aureus_subsp._aureus_NCTC_8325_GCF_000013425.1.fna.gz',
  'Stenotrophomonas_maltophilia_GCF_900475405.1.fna.gz',
  'Streptococcus_anginosus_GCF_001412635.1.fna.gz',
  'Streptococcus_dysgalactiae_GCF_016724885.1.fna.gz',
  'Streptococcus_pneumoniae_GCF_002076835.1.fna.gz',
  'Streptococcus_pseudopneumoniae_IS7493_GCF_000221985.1.fna.gz',
  'Streptococcus_pyogenes_GCF_900475035.1.fna.gz',
  'Vibrio_alginolyticus_GCA_023650915.1.fna.gz',
  'Vibrio_alginolyticus_GCF_009665435.1.fna.gz',
  'Vibrio_cholerae_GCA_009665515.2.fna.gz',
  'Vibrio_cholerae_GCF_008369605.1.fna.gz',
  'Vibrio_cidicii_GCA_009665415.1.fna.gz',
  'Vibrio_cincinnatiensis_GCF_009665395.1.fna.gz',
  'Vibrio_fluvialis_GCF_009665355.1.fna.gz',
  'Vibrio_furnissii_GCF_009665335.1.fna.gz',
  'Vibrio_harveyi_GCF_009665315.1.fna.gz',
  'Vibrio_metoecus_GCF_009665255.1.fna.gz',
  'Vibrio_metoecus_GCF_009665275.1.fna.gz',
  'Vibrio_metschnikovii_GCF_009665235.1.fna.gz',
  'Vibrio_mimicus_GCF_009665195.1.fna.gz',
  'Vibrio_mimicus_MB451_GCF_000176375.1.fna.gz',
  'Vibrio_navarrensis_GCF_009665215.1.fna.gz',
  'Vibrio_navarrensis_GCF_012275065.1.fna.gz',
  'Vibrio_paracholerae_GCA_003311965.1.fna.gz',
  'Vibrio_parahaemolyticus_GCF_009665495.1.fna.gz',
  'Vibrio_parahaemolyticus_RIMD_2210633_GCF_000196095.1.fna.gz',
  'Vibrio_vulnificus_GCF_002204915.1.fna.gz',
  'Vibrio_vulnificus_GCF_009665455.1.fna.gz',
  'Vibrio_vulnificus_GCF_009665475.1.fna.gz'
]

// Getting the reference genomes for fastANI
ch_fastani_genomes = Channel.empty()

for (genome in included_genomes) {
  ch_fastani_genomes = ch_fastani_genomes.mix(Channel.fromPath("${workflow.projectDir}/db/" + genome, type: "file"))
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for input files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// using a sample sheet with the column header pf 'sample,fastq_1,fastq_2'
ch_input_reads = params.sample_sheet
  ? Channel
    .fromPath("${params.sample_sheet}", type: "file")
    .view { "Sample sheet found : ${it}" }
    .splitCsv( header: true, sep: ',' )
    .map { row -> tuple( "${row.sample}", [ file("${row.fastq_1}"), file("${row.fastq_2}") ]) }
  : Channel.empty()

// using a sample sheet for fasta files (no header)
ch_input_fastas = params.fasta_list
  ? Channel
    .fromPath("${params.fasta_list}", type: "file")
    .view { "Fasta list found : ${it}" }
    .splitText()
    .map{ it -> it.trim()}
    .map{ it -> file(it) }
    .map{ it -> tuple(it.baseName, it) }
  : Channel.empty()

// Getting the fastq files
Channel
  .fromFilePairs(["${params.reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}",
                  "${params.reads}/*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}"], size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .view { "Paired-end fastq files found : ${it[0]}" }
  .unique()
  .set { ch_reads }

// Getting contig or fasta files
Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
  .map { file -> tuple(file.baseName, file) }
  .unique()
  .set { ch_fastas }

// Getting fasta files that have been annotated with prokka
Channel
  .fromPath("${params.gff}/*.gff", type: "file")
  .view { "gff file : $it" }
  .unique()
  .set { ch_gffs }

// Getting accession for downloading
ch_sra_accessions = Channel.from( params.sra_accessions )

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Channels for database files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Getting the file with genome sizes of common organisms for fastqcscan. The End User can use their own file and set with a param
Channel
  .fromPath(params.genome_sizes, type: "file")
  .ifEmpty{
    println("The genome sizes file for this workflow are missing!")
    exit 1}
  .set { ch_genome_sizes }

// Getting the database for blobtools
ch_blast_db = params.blast_db
  ? Channel
    .fromPath(params.blast_db, type: "dir")
    .ifEmpty{
      println("No blast database was found at ${params.blast_db}")
      println("Set 'params.blast_db' to directory with blast database")
      exit 1
      }
      .view { "Local Blast Database for Blobtools : $it" }
    : Channel.empty()

// Getting the kraken2 database
ch_kraken2_db = params.kraken2_db
  ? Channel
    .fromPath(params.kraken2_db, type: "dir")
    .ifEmpty{
      println("No kraken2 database was found at ${params.kraken2_db}")
      println("Set 'params.kraken2_db' to directory with kraken2 database")
      exit 1
      }
      .view { "Local kraken2 database : $it" }
  : Channel.empty()

// Getting the mash reference
ch_mash_db = params.mash_db 
  ? Channel
    .fromPath(params.mash_db, type: "file")
    .ifEmpty{
      println("No mash database was found at ${params.mash_db}")
      println("Set 'params.mash_db' to file of pre-sketched mash reference")
      exit 1
      }
    .view { "Mash reference : $it" }
  : Channel.empty()

//# user supplied fastani reference genomes
if ( params.fastani_ref ) {
  Channel
    .of( params.fastani_ref )
    .splitCsv()
    .flatten()
    .map { it -> file(it) }
    .view{ "Additional fastani reference genomes : $it" }
    .set { ch_fastani_genomes_input }

  ch_fastani_genomes = ch_fastani_genomes.mix(ch_fastani_genomes_input)
}

if ( params.fastani_ref_list ) {
  Channel.fromPath(params.fastani_ref_list, type: "file")
    .splitText()
    .map( it -> it.trim())
    .map{ it -> file(it) }
    .view{ "Additional fastani reference genome from file : $it" }
    .set{ ch_fastani_ref_list }

  ch_fastani_genomes = ch_fastani_genomes.mix(ch_fastani_ref_list)
}

println("The files and directory for results is " + params.outdir )
println("The maximum number of CPUS for any one process is ${params.maxcpus}")

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Workflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow {

  ch_for_multiqc   = Channel.empty()
  ch_for_summary   = Channel.empty()
  ch_for_flag      = Channel.empty()
  ch_top_hit       = Channel.empty()

  // getting test files
  if ( ! params.sra_accessions.isEmpty() ) { 
    test(ch_sra_accessions)
    ch_raw_reads = ch_reads.mix(test.out.fastq).mix(ch_input_reads)
  } else {
    ch_raw_reads = ch_reads.mix(ch_input_reads)
  }

  // running DONUT FALLS first (under development)
  // if ( params.donut_falls_wf )  { donut_falls(input) }

  // either phoenix or de_novo_alignment is required
  de_novo_alignment(ch_raw_reads)
  
  ch_contigs       = ch_fastas.mix(de_novo_alignment.out.contigs).mix(ch_input_fastas)
  ch_clean_reads   = de_novo_alignment.out.clean_reads

  ch_for_multiqc   = ch_for_multiqc.mix(de_novo_alignment.out.for_multiqc)

  // optional subworkflow blobtools (useful for interspecies contamination)
  if ( params.blast_db ) {
    blobtools(ch_clean_reads, ch_contigs, ch_blast_db )

    ch_for_multiqc = ch_for_multiqc.mix(blobtools.out.for_multiqc)
    ch_for_summary = ch_for_summary.mix(blobtools.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(blobtools.out.for_flag)
  }

  // optional subworkflow kraken2 (useful for interspecies contamination)
  if ( params.kraken2_db ) {
    kmer_taxonomic_classification(ch_clean_reads, ch_kraken2_db )

    ch_for_multiqc = ch_for_multiqc.mix(kmer_taxonomic_classification.out.for_multiqc)
    ch_for_summary = ch_for_summary.mix(kmer_taxonomic_classification.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(kmer_taxonomic_classification.out.for_flag)
  } 
  
  if (params.extras ) {
    // subworkflow mash for species determination
    min_hash_distance(ch_clean_reads, ch_contigs, ch_mash_db)

    ch_for_summary = ch_for_summary.mix(min_hash_distance.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(min_hash_distance.out.for_flag)

    // determining organisms in sample
    average_nucleotide_identity(
      ch_for_summary.collect(),
      ch_contigs,
      ch_fastani_genomes,
      dataset_script)

    ch_for_summary = ch_for_summary.mix(average_nucleotide_identity.out.for_summary)
    ch_for_flag    = ch_for_flag.mix(average_nucleotide_identity.out.for_flag)
    ch_top_hit     = ch_top_hit.mix(average_nucleotide_identity.out.top_hit)
    ch_datasets    = average_nucleotide_identity.out.datasets_summary.ifEmpty('none')

    ch_contigs
      .join(min_hash_distance.out.mash_err)
      .join(min_hash_distance.out.for_flag)
      .join(average_nucleotide_identity.out.for_flag)
      .combine(ch_genome_sizes)
      .combine(ch_datasets)
      .set{ ch_size }

    // getting all the other information
    information(
      ch_raw_reads, 
      ch_contigs, 
      ch_for_flag, 
      ch_size,
      summfle_script)

    ch_for_summary = ch_for_summary.mix(information.out.for_summary)
    ch_for_multiqc = ch_for_multiqc.mix(information.out.for_multiqc)
  } 

  // optional subworkflow for comparing shared genes
  if ( params.msa ) {
    phylogenetic_analysis(
      evaluat_script, 
      ch_contigs.ifEmpty([]), 
      ch_gffs.ifEmpty([]), 
      ch_top_hit.ifEmpty([]))
    
    ch_for_multiqc   = ch_for_multiqc.mix(phylogenetic_analysis.out.for_multiqc)
  }

  // getting a summary of everything
  if (params.extras ) { 
    report(
      ch_raw_reads, 
      ch_fastas, 
      ch_for_multiqc.collect(), 
      ch_for_summary.concat(summary_script).collect())
  }
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Final Steps

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("MultiQC report can be found at ${params.outdir}/multiqc/multiqc_report.html")
    println("Summary can be found at ${params.outdir}/grandeur_summary.tsv")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}

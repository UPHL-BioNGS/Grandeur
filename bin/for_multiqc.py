#!/usr/bin/env python

import shutil
from os.path import exists
import pandas as pd
import numpy as np

##########################################
# defining files                         #
##########################################

# input files
amrfinder_input         = 'amrfinderplus.txt'
amrfinder_output        = 'amrfinderplus_mqc.txt'
blobtools_input         = 'blobtools_summary.txt'
blobtools_output        = 'blobtools_mqc.tsv'
circulocov_input        = 'circulocov_summary.tsv'
circulocov_output       = 'circulocov_mqc.tsv'
core_genome_input       = 'core_genome_evaluation.csv'
core_genome_output      = 'core_genome_evaluation_mqc.csv' 
drprg_input             = 'drprg_summary.tsv'
drprg_output            = 'drprg_mqc.tsv'
elgato_input            = 'elgato_summary.tsv'
elgato_output           = 'elgato_mqc.tsv'
emmtyper_input          = 'emmtyper_summary.tsv'
emmtyper_output         = 'emmtyper_mqc.tsv'
fastani_input           = 'fastani_summary.csv'
fastani_output          = 'fastani_mqc.csv'
heatcluster_input       = 'heatcluster.png'
heatcluster_output      = 'heatcluster_mqc.png'
kaptive_input           = 'kaptive_summary.txt'
kaptive_output          = 'kaptive_mqc.tsv'
kleborate_input         = 'kleborate_results.tsv'
kleborate_output        = 'kleborate_mqc.tsv'
mash_input              = 'mash_summary.csv'
mash_output             = 'mash_mqc.csv'
mlst_input              = 'mlst_summary.tsv'
mlst_output             = 'mlst_mqc.tsv'
mykrobe_input           = 'mykrobe_summary.csv'
mykrobe_output          = 'mykrobe_summary_mqc.csv'
pbp_input               = 'pbptyper_summary.tsv'
pbp_output              = 'pbptyper_mqc.tsv'
phytreeviz_iqtr_input   = 'iqtree_tree.png'
phytreeviz_iqtr_output  = 'phytreeviz_iqtree2_mqc.png' 
phytreeviz_mshtr_input  = 'mashtree_tree.png'
phytreeviz_mshtr_output = 'phytreeviz_mashtree_mqc.png'
plasmidfinder_input     = 'plasmidfinder_result.tsv'
plasmidfinder_output    = 'plasmidfinder_mqc.tsv' 
seqsero2_input          = 'seqsero2_results.txt'
seqsero2_output         = 'seqsero2_mqc.txt'
serotypefinder_input    = 'serotypefinder_results.txt'
serotypefinder_output   = 'serotypefinder_mqc.txt'
shigatyper_input        = 'shigatyper_summary.txt'
shigatyper_output       = 'shigatyper_mqc.txt'
shigatyper_hit_input    = 'shigatyper_hits.txt' 
shigatyper_hit_output   = 'shigatyper_hits_mqc.tsv'
snpdists_input          = 'snp_matrix.txt'
snpdists_output         = 'snpdists_matrix_mqc.txt' 

##########################################
# getting ready for multiqc              #
##########################################

if exists(blobtools_input) :
    line = "# plot_type: 'bargraph'\n"
    line = line + "# section_name: 'Blobtools'\n"
    line = line + "# description: 'Visualisation, quality control and taxonomic partitioning of genome datasets'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    with open(blobtools_output,'w') as f:
        f.write(line)

    blobtools_df = pd.read_table(blobtools_input)
    blobtools_df = blobtools_df[~blobtools_df['name'].isin(['all', 'no-hit', 'undef'])]

    samples = blobtools_df['sample'].drop_duplicates().tolist()
    organisms = sorted(blobtools_df['name'].drop_duplicates().tolist())
    blobtools_result_df = pd.DataFrame(columns=["sample"] + organisms)

    for sample in samples:
        result = [sample]
        for organism in organisms:
            readper = blobtools_df.loc[(blobtools_df['sample'] == sample) & (blobtools_df['name'] == organism), 'bam0_read_map_p']
            orgper = readper.iloc[0] if not readper.empty else 0
            result = result + [orgper]
        
        blobtools_result_df.loc[len(blobtools_result_df.index)] = result

    blobtools_result_df.to_csv(blobtools_output, index=False, sep="\t", mode='a')

if exists(drprg_input):
    # sample  genes_present   genes_absent    susceptibility
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'Dr. PRG'\n"
    line = line + "# description: 'Drug resistance Prediction with Reference Graphs'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# format: 'tsv'\n"
    with open(drprg_output,'w') as f:
        f.write(line)

    drprg_df = pd.read_table(drprg_input)
    drprg_df.to_csv(drprg_output, index=False, sep="\t", mode='a')

if exists(circulocov_input):
    # sample  circ    contigs length  illumina_numreads       illumina_covbases       illumina_coverage       illumina_meandepth
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'CirculoCov'\n"
    line = line + "# description: 'Coverage Estimates through Mapping'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     length:\n"
    line = line + "#         title: 'length'\n"
    line = line + "#         description: 'Total length'\n"
    line = line + "#     illumina_numreads:\n"
    line = line + "#         title: 'numreads'\n"
    line = line + "#         description: 'illumina_numreads'\n"
    line = line + "#     illumina_covbases:\n"
    line = line + "#         title: 'covbases'\n"
    line = line + "#         description: 'illumina_covbases'\n"
    line = line + "#     illumina_coverage:\n"
    line = line + "#         title: 'coverage'\n"
    line = line + "#         description: 'illumina_coverage'\n"
    line = line + "#     illumina_meandepth:\n"
    line = line + "#         title: 'meandepth'\n"
    line = line + "#         description: 'illumina_meandepth'\n"
    with open(circulocov_output,'w') as f:
        f.write(line)

    circulocov_df = pd.read_table(circulocov_input)
    circulocov_df = circulocov_df[circulocov_df['contigs'] == 'all']
    circulocov_df = circulocov_df.drop(['circ', 'contigs'], axis=1)
    circulocov_df.to_csv(circulocov_output, index=False, sep="\t", mode='a')

if exists(elgato_input):
    # Sample  ST      flaA    pilE    asd     mip     mompS   proA    neuA_neuAH
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'el_gato'\n"
    line = line + "# description: 'Epidemiology of Legionella : Genome-bAsed Typing'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     ST:\n"
    line = line + "#         title: 'ST'\n"
    line = line + "#         description: 'ST'\n"
    with open(elgato_output,'w') as f:
        f.write(line)

    elgato_df = pd.read_table(elgato_input)
    elgato_df.to_csv(elgato_output, index=False, sep="\t", mode='a')

if exists(emmtyper_input):
    # sample  Isolate name    Number of BLAST hits    Number of clusters      Predicted emm-type      Position(s)     Possible emm-like alleles       emm-like position(s)    EMM cluster
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'emmtyper'\n"
    line = line + "# description: 'Emm Automatic Isolate Labeller'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     Number of BLAST hits:\n"
    line = line + "#         title: 'Number of BLAST hits'\n"
    line = line + "#         description: 'Number of BLAST hits'\n"
    line = line + "#     Number of clusters:\n"
    line = line + "#         title: 'Number of clusters'\n"
    line = line + "#         description: 'Number of clusters'\n"
    line = line + "#     Predicted emm-type:\n"
    line = line + "#         title: 'Predicted emm-type'\n"
    line = line + "#         description: 'Predicted emm-type'\n"
    line = line + "#     Position(s):\n"
    line = line + "#         title: 'Position(s)'\n"
    line = line + "#         description: 'Position(s)'\n"
    line = line + "#     Possible emm-like alleles:\n"
    line = line + "#         title: 'Possible emm-like alleles'\n"
    line = line + "#         description: 'Possible emm-like alleles'\n"
    line = line + "#     emm-like position(s):\n"
    line = line + "#         title: 'emm-like position(s)'\n"
    line = line + "#         description: 'emm-like position(s)'\n"
    line = line + "#     EMM cluster:\n"
    line = line + "#         title: 'EMM cluster'\n"
    line = line + "#         description: 'EMM cluster'\n"
    with open(emmtyper_output,'w') as f:
        f.write(line)

    emmtyper_df = pd.read_table(emmtyper_input)
    emmtyper_df.to_csv(emmtyper_output, index=False, sep="\t", mode='a')

if exists(mash_input) :
    line = "# plot_type: 'bargraph'\n"
    line = line + "# section_name: 'Mash dist'\n"
    line = line + "# description: 'Fast genome and metagenome distance estimation using MinHash'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    with open(mash_output,'w') as f:
        f.write(line)

    mash_df = pd.read_csv(mash_input)
    samples = mash_df['sample'].drop_duplicates().tolist()
    organisms = sorted(mash_df['organism'].drop_duplicates().tolist())
    mash_result_df = pd.DataFrame(columns=["sample"] + organisms)

    for sample in samples:
        df_len = len(mash_result_df)
        mash_result_df.loc[df_len] = pd.Series(dtype='str')
        mash_result_df.at[df_len, "sample"] = sample
        sample_df = mash_df[mash_df['sample'] == sample].copy()

        for organism in organisms:
            mash_result_df.at[df_len, organism] = (sample_df['organism'] == organism).sum()
        #mash_result_df = mash_result_df.reset_index()

        # only worked in python 3.12 :(
        #sample_df = mash_df[mash_df['sample' ] == sample].copy()
        #counts_df = sample_df['organism'].value_counts().reset_index()
        #counts_df = counts_df.rename(columns={"index": "organism", 0: "count"}) 
        #counts_df = counts_df.set_index('organism')
        #counts_df.columns = sample
        #counts_df = counts_df.transpose()
        #mash_result_df = pd.concat([mash_result_df, counts_df], axis=0, join='outer')

    mash_result_df = mash_result_df.astype(str).fillna("0")
    mash_result_df.to_csv(mash_output, index=False,mode='a')

if exists(fastani_input):
    # sample,query,reference,ANI estimate,total query sequence fragments,fragments aligned as orthologous matches
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'FastANI'\n"
    line = line + "# description: 'Fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     sample:\n"
    line = line + "#         title: 'sample'\n"
    line = line + "#         description: 'sample'\n"
    line = line + "#     reference:\n"
    line = line + "#         title: 'reference'\n"
    line = line + "#         description: 'reference'\n"
    line = line + "#     ANI estimate:\n"
    line = line + "#         title: 'ANI estimate'\n"
    line = line + "#         description: 'ANI estimate'\n"
    line = line + "#     total query sequence:\n"
    line = line + "#         title: 'total query sequence'\n"
    line = line + "#         description: 'total query sequence'\n"
    line = line + "#     fragments aligned as orthologous matches:\n"
    line = line + "#         title: 'fragments aligned as orthologous matches'\n"
    line = line + "#         description: 'fragments aligned as orthologous matches'\n"
    with open(fastani_output,'w') as f:
        f.write(line)

    fastani_df = pd.read_csv(fastani_input)
    fastani_df = fastani_df[fastani_df['ANI estimate'] > 90]
    fastani_df.to_csv(fastani_output, mode='a')

if exists(shigatyper_hit_input):
    # sample          Hit     Number of reads Length Covered  reference length        % covered       Number of variants      % accuracy
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'ShigaTyper Hits'\n"
    line = line + "# description: 'ShigaTyper Hits'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     sample:\n"
    line = line + "#         title: 'sample'\n"
    line = line + "#         description: 'sample'\n"
    line = line + "#     Hit:\n"
    line = line + "#         title: 'Hit'\n"
    line = line + "#         description: 'Hit'\n"
    with open(shigatyper_hit_output,'w') as f:
        f.write(line)

    shigatyper_df = pd.read_table(shigatyper_hit_input)
    shigatyper_df = shigatyper_df.drop(shigatyper_df.columns[1], axis=1)
    shigatyper_df = shigatyper_df.iloc[:, :2]
    shigatyper_df.to_csv(shigatyper_hit_output, sep="\t", mode='a')
    #if [ -f 'shigatyper_results.txt' ]     ; then awk '{print \$1 "_" \$2 "\t" \$3}' shigatyper_results.txt > shigatyper_mqc.tsv ; fi

if exists(shigatyper_input):
    # sample  prediction      ipaB    notes
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'ShigaTyper'\n"
    line = line + "# description: 'Shigella serotype'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     prediction:\n"
    line = line + "#         title: 'prediction'\n"
    line = line + "#         description: 'prediction'\n"
    line = line + "#     ipaB:\n"
    line = line + "#         title: 'ipaB'\n"
    line = line + "#         description: 'ipaB'\n"
    line = line + "#     notes:\n"
    line = line + "#         title: 'notes'\n"
    line = line + "#         description: 'notes'\n"
    with open(shigatyper_output,'w') as f:
        f.write(line)

    shigatyper_df = pd.read_table(shigatyper_input)
    shigatyper_df.to_csv(shigatyper_output, index=False, sep="\t", mode='a')

if exists(amrfinder_input):
    # Name    Protein identifier      Contig id       Start   Stop    Strand  Gene symbol     Sequence name   Scope   Element type    Element subtype Class   Subclass        Method  Target length        Reference sequence length       % Coverage of reference sequence        % Identity to reference sequence        Alignment length        Accession of closest sequence   Name of closest sequence     HMM id  HMM description
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'AMRFinderPlus'\n"
    line = line + "# description: 'NCBI Antimicrobial Resistance Gene Finder'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     Name:\n"
    line = line + "#         title: 'Name'\n"
    line = line + "#         description: 'Name'\n"
    line = line + "#     Contig id:\n"
    line = line + "#         title: 'Contig id'\n"
    line = line + "#         description: 'Contig gene was found on'\n"
    line = line + "#     Start:\n"
    line = line + "#         title: 'Start'\n"
    line = line + "#         description: 'Start position of AMR gene'\n"
    line = line + "#     Stop:\n"
    line = line + "#         title: 'Stop'\n"
    line = line + "#         description: 'Stop position of AMR gene'\n"
    # line = line + "#     Strand:\n"
    # line = line + "#         title: 'Strand'\n"
    # line = line + "#         description: 'Strand of AMR gene'\n"
    line = line + "#     Gene symbol:\n"
    line = line + "#         title: 'Gene symbol'\n"
    line = line + "#         description: 'Gene symbol'\n"
    # line = line + "#     Sequence name:\n"
    # line = line + "#         title: 'Sequence name'\n"
    # line = line + "#         description: 'Sequence name'\n"
    line = line + "#     Scope:\n"
    line = line + "#         title: 'Scope'\n"
    line = line + "#         description: 'Scope'\n"
    line = line + "#     Element type:\n"
    line = line + "#         title: 'Element type'\n"
    line = line + "#         description: 'Element type'\n"
    line = line + "#     Element subtype:\n"
    line = line + "#         title: 'Element subtype'\n"
    line = line + "#         description: 'Element subtype'\n"
    line = line + "#     Class:\n"
    line = line + "#         title: 'Class'\n"
    line = line + "#         description: 'Class'\n"
    line = line + "#     Subclass:\n"
    line = line + "#         title: 'Subclass'\n"
    line = line + "#         description: 'Subclass'\n"
    # line = line + "#     Method:\n"
    # line = line + "#         title: 'Method'\n"
    # line = line + "#         description: 'Method'\n"
    # line = line + "#     Target length:\n"
    # line = line + "#         title: 'Target length'\n"
    # line = line + "#         description: 'Target length'\n"
    # line = line + "#     Reference sequence length:\n"
    # line = line + "#         title: 'Reference sequence length'\n"
    # line = line + "#         description: 'Reference sequence length'\n"
    line = line + "#     Per Coverage of reference sequence:\n"
    line = line + "#         title: 'Per Coverage'\n"
    line = line + "#         description: 'Percent Coverage of reference sequence'\n"
    line = line + "#     Per Identity:\n"
    line = line + "#         title: 'Per Identity to reference sequence'\n"
    line = line + "#         description: 'Percent Identity to reference sequence'\n"
    # line = line + "#     Alignment length:\n"
    # line = line + "#         title: 'Alignment length'\n"
    # line = line + "#         description: 'Alignment length'\n"
    line = line + "#     Accession of closest sequence:\n"
    line = line + "#         title: 'Accession of closest sequence'\n"
    line = line + "#         description: 'Accession of closest sequence'\n"
    line = line + "#     Name of closest sequence:\n"
    line = line + "#         title: 'Name of closest sequence'\n"
    line = line + "#         description: 'Name of closest sequence'\n"
    with open(amrfinder_output,'w') as f:
        f.write(line)

    amrfinder_df = pd.read_table(amrfinder_input)
    amrfinder_df = amrfinder_df.replace('%', 'Per', regex=True)
    #amrfinder_df = amrfinder_df.replace(' ', '_', regex=True)
    amrfinder_df.to_csv(amrfinder_output, sep="\t", mode='a')

if exists(kaptive_input):
    # Assembly        Best match locus        Best match type Match confidence        Problems        Coverage        Identity        Length discrepancy      Expected genes in locus Expected genes in locus, details        Missing expected genes       Other genes in locus    Other genes in locus, details   Expected genes outside locus    Expected genes outside locus, details   Other genes outside locus       Other genes outside locus, details
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'Kaptive'\n"
    line = line + "# description: 'Vibrio species complex'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     Assembly:\n"
    line = line + "#         title: 'Assembly'\n"
    line = line + "#         description: 'Assembly'\n"
    line = line + "#     Best match locus:\n"
    line = line + "#         title: 'Best match locus'\n"
    line = line + "#         description: 'Best match locus'\n"
    line = line + "#     Best match type:\n"
    line = line + "#         title: 'Best match type'\n"
    line = line + "#         description: 'Best match type'\n"
    line = line + "#     Match confidence:\n"
    line = line + "#         title: 'Match confidence'\n"
    line = line + "#         description: 'Match confidence'\n"
    line = line + "#     Problems:\n"
    line = line + "#         title: 'Problems'\n"
    line = line + "#         description: 'Problems'\n"
    line = line + "#     Coverage:\n"
    line = line + "#         title: 'Coverage'\n"
    line = line + "#         description: 'Coverage'\n"
    line = line + "#     Identity:\n"
    line = line + "#         title: 'Identity'\n"
    line = line + "#         description: 'Identity'\n"
    line = line + "#     Length discrepancy:\n"
    line = line + "#         title: 'Length discrepancy'\n"
    line = line + "#         description: 'Length discrepancy'\n"
    line = line + "#     Expected genes in locus:\n"
    line = line + "#         title: 'Expected genes in locus'\n"
    line = line + "#         description: 'Expected genes in locus'\n"
    line = line + "#     Expected genes in locus, details:\n"
    line = line + "#         title: 'Expected genes in locus, details'\n"
    line = line + "#         description: 'Expected genes in locus, details'\n"
    line = line + "#     Missing expected genes:\n"
    line = line + "#         title: 'Missing expected genes'\n"
    line = line + "#         description: 'Missing expected genes'\n"
    line = line + "#     Other genes in locus:\n"
    line = line + "#         title: 'Other genes in locus'\n"
    line = line + "#         description: 'Other genes in locus'\n"
    line = line + "#     Other genes in locus, details:\n"
    line = line + "#         title: 'Other genes in locus, details'\n"
    line = line + "#         description: 'Other genes in locus, details'\n"
    line = line + "#     Expected genes outside locus:\n"
    line = line + "#         title: 'Expected genes outside locus'\n"
    line = line + "#         description: 'Expected genes outside locus'\n"
    line = line + "#     Expected genes outside locus, details:\n"
    line = line + "#         title: 'Expected genes outside locus, details'\n"
    line = line + "#         description: 'Expected genes outside locus, details'\n"
    line = line + "#     Other genes outside locus:\n"
    line = line + "#         title: 'Other genes outside locus'\n"
    line = line + "#         description: 'Other genes outside locus'\n"
    line = line + "#     Other genes outside locus, details:\n"
    line = line + "#         title: 'Other genes outside locus, details'\n"
    line = line + "#         description: 'Other genes outside locus, details'\n"
    with open(kaptive_output,'w') as f:
        f.write(line)

    kaptive_df = pd.read_table(kaptive_input)
    kaptive_df = kaptive_df.replace(' ', '_', regex=True)
    kaptive_df.to_csv(kaptive_output, index=True, sep="\t", mode='a')


if exists(kleborate_input):
    # sample  strain  species species_match   contig_count    N50     largest_contig  total_size      ambiguous_bases QC_warnings     ST      virulence_score resistance_score        num_resistance_classes  num_resistance_genes    Yersiniabactin       YbST    Colibactin      CbST    Aerobactin      AbST    Salmochelin     SmST    RmpADC  RmST    rmpA2   wzi     K_locus K_type  K_locus_problems        K_locus_confidence      K_locus_identity        K_locus_missing_genes        O_locus O_type  O_locus_problems        O_locus_confidence      O_locus_identity        O_locus_missing_genes   AGly_acquired   Col_acquired    Fcyn_acquired   Flq_acquired    Gly_acquired    MLS_acquired    Phe_acquired Rif_acquired    Sul_acquired    Tet_acquired    Tgc_acquired    Tmt_acquired    Bla_acquired    Bla_inhR_acquired       Bla_ESBL_acquired       Bla_ESBL_inhR_acquired  Bla_Carb_acquired       Bla_chr SHV_mutations   Omp_mutations        Col_mutations   Flq_mutations   truncated_resistance_hits       spurious_resistance_hits        Chr_ST  gapA    infB    mdh     pgi     phoE    rpoB    tonB    ybtS    ybtX    ybtQ    ybtP    ybtA    irp2    irp1ybtU     ybtT    ybtE    fyuA    clbA    clbB    clbC    clbD    clbE    clbF    clbG    clbH    clbI    clbL    clbM    clbN    clbO    clbP    clbQ    iucA    iucB    iucC    iucD    iutA    iroB    iroC    iroD    iroN    rmpArmpD     rmpC    spurious_virulence_hits
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'Kleborate'\n"
    line = line + "# description: 'Virulence Scoring'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     strain:\n"
    line = line + "#         title: 'strain'\n"
    line = line + "#         description: 'strain'\n"
    line = line + "#     species:\n"
    line = line + "#         title: 'species'\n"
    line = line + "#         description: 'species'\n"
    line = line + "#     species_match:\n"
    line = line + "#         title: 'species_match'\n"
    line = line + "#         description: 'species_match'\n"
    line = line + "#     virulence_score:\n"
    line = line + "#         title: 'virulence_score'\n"
    line = line + "#         description: 'virulence_score'\n"
    line = line + "#     resistance_score:\n"
    line = line + "#         title: 'resistance_score'\n"
    line = line + "#         description: 'resistance_score'\n"
    line = line + "#     contig_count:\n"
    line = line + "#         title: 'contig_count'\n"
    line = line + "#         description: 'contig_count'\n"
    line = line + "#     N50:\n"
    line = line + "#         title: 'N50'\n"
    line = line + "#         description: 'N50'\n"
    line = line + "#     largest_contig:\n"
    line = line + "#         title: 'largest_contig'\n"
    line = line + "#         description: 'largest_contig'\n"
    line = line + "#     total_size:\n"
    line = line + "#         title: 'total_size'\n"
    line = line + "#         description: 'total_size'\n"
    line = line + "#     ambiguous_bases:\n"
    line = line + "#         title: 'ambiguous_bases'\n"
    line = line + "#         description: 'ambiguous_bases'\n"
    line = line + "#     QC_warnings:\n"
    line = line + "#         title: 'QC_warnings'\n"
    line = line + "#         description: 'QC_warnings'\n"
    line = line + "#     ST:\n"
    line = line + "#         title: 'ST'\n"
    line = line + "#         description: 'ST'\n"
    with open(kleborate_output,'w') as f:
        f.write(line)

    kleborate_df = pd.read_table(kleborate_input)
    kleborate_df = kleborate_df.iloc[:, [1] + list(range(2, 12))]
    kleborate_df.to_csv(kleborate_output, index=False, sep="\t", mode='a')
    #if [ -f 'kleborate_results.tsv' ]      ; then cut -f 1,3-12 kleborate_results.tsv > kleborate_mqc.tsv       ; fi

if exists(mlst_input):
    # sample  filename        matching PubMLST scheme ST      ID1     ID2     ID3     ID4     ID5     ID6     ID7     ID8     ID9     ID10    ID11    ID12    ID13    ID14    ID15
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'mlst'\n"
    line = line + "# description: 'Scan contig files against traditional PubMLST typing schemes'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     matching PubMLST scheme:\n"
    line = line + "#         title: 'matching PubMLST scheme'\n"
    line = line + "#         description: 'matching PubMLST scheme'\n"
    line = line + "#     ST:\n"
    line = line + "#         title: 'ST'\n"
    line = line + "#         description: 'ST'\n"
    with open(mlst_output,'w') as f:
        f.write(line)

    mlst_df = pd.read_table(mlst_input)
    mlst_df = mlst_df.replace(' ', '_', regex=True)
    mlst_df = mlst_df.loc[:, ['sample', 'matching PubMLST scheme', 'ST']] 
    mlst_df.to_csv(mlst_output, index=False, sep="\t", mode='a')

if exists(mykrobe_input):
    # "sample","drug","susceptibility","variants (dna_variant-AA_variant:ref_kmer_count:alt_kmer_count:conf) [use --format json for more info]","genes (prot_mut-ref_mut:percent_covg:depth) [use --format json for more info]","mykrobe_version","files","probe_sets","genotype_model","kmer_size","phylo_group","species","lineage","phylo_group_per_covg","species_per_covg","lineage_per_covg","phylo_group_depth","species_depth","lineage_depth"
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'Mykrobe'\n"
    line = line + "# description: 'Mycobacteria taxonomy and names'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     drug:\n"
    line = line + "#         title: 'drug'\n"
    line = line + "#         description: 'drug'\n"
    line = line + "#     susceptibility:\n"
    line = line + "#         title: 'susceptibility'\n"
    line = line + "#         description: 'susceptibility'\n"
    line = line + "#     phylo_group:\n"
    line = line + "#         title: 'phylo_group'\n"
    line = line + "#         description: 'phylo_group'\n"
    line = line + "#     species:\n"
    line = line + "#         title: 'species'\n"
    line = line + "#         description: 'species'\n"
    line = line + "#     lineage:\n"
    line = line + "#         title: 'lineage'\n"
    line = line + "#         description: 'lineage'\n"
    with open(mykrobe_output,'w') as f:
        f.write(line)

    mykrobe_df = pd.read_csv(mykrobe_input)
    mykrobe_df.to_csv(mykrobe_output, index=False, mode='a')

if exists (pbp_input):
    # sample  pbptype ani     1A_coverage     1A_pident       1A_bitscore     2B_coverage     2B_pident       2B_bitscore     2X_coverage     2X_pident       2X_bitscore     comment
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'pbptyper'\n"
    line = line + "# description: 'In silico Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae assemblies'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     pbptype:\n"
    line = line + "#         title: 'pbptype'\n"
    line = line + "#         description: 'pbptype'\n"
    line = line + "#     ani:\n"
    line = line + "#         title: 'ani'\n"
    line = line + "#         description: 'ani'\n"
    line = line + "#     1A_coverage:\n"
    line = line + "#         title: '1A_coverage'\n"
    line = line + "#         description: '1A_coverage'\n"
    line = line + "#     1A_pident:\n"
    line = line + "#         title: '1A_pident'\n"
    line = line + "#         description: '1A_pident'\n"
    line = line + "#     1A_bitscore:\n"
    line = line + "#         title: '1A_bitscore'\n"
    line = line + "#         description: '1A_bitscore'\n"
    line = line + "#     2B_coverage:\n"
    line = line + "#         title: '2B_coverage'\n"
    line = line + "#         description: '2B_coverage'\n"
    line = line + "#     2B_pident:\n"
    line = line + "#         title: '2B_pident'\n"
    line = line + "#         description: '2B_pident'\n"
    line = line + "#     2B_bitscore:\n"
    line = line + "#         title: '2B_bitscore'\n"
    line = line + "#         description: '2B_bitscore'\n"
    line = line + "#     2X_coverage:\n"
    line = line + "#         title: '2X_coverage'\n"
    line = line + "#         description: '2X_coverage'\n"
    line = line + "#     2X_pident:\n"
    line = line + "#         title: '2X_pident'\n"
    line = line + "#         description: '2X_pident'\n"
    line = line + "#     2X_bitscore:\n"
    line = line + "#         title: '2X_bitscore'\n"
    line = line + "#         description: '2X_bitscore'\n"
    line = line + "#     comment:\n"
    line = line + "#         title: 'comment'\n"
    line = line + "#         description: 'comment'\n"
    with open(pbp_output,'w') as f:
        f.write(line)

    pbp_df = pd.read_table(pbp_input)
    pbp_df.to_csv(pbp_output, index=False, sep="\t", mode='a')

if exists(plasmidfinder_input):
    # sample  Database        Plasmid Identity        Query / Template length Contig  Position in contig      Note    Accession number
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'PlasmidFinder'\n"
    line = line + "# description: 'Identifies plasmids in total or partial sequenced isolates of bacteria.n'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     Database:\n"
    line = line + "#         title: 'Database'\n"
    line = line + "#         description: 'Database'\n"
    line = line + "#     Plasmid Identity:\n"
    line = line + "#         title: 'Plasmid Identity'\n"
    line = line + "#         description: 'Plasmid Identity'\n"
    line = line + "#     Query / Template length:\n"
    line = line + "#         title: 'Query / Template length'\n"
    line = line + "#         description: 'Query / Template length'\n"
    line = line + "#     Contig:\n"
    line = line + "#         title: 'Contig'\n"
    line = line + "#         description: 'Contig'\n"
    line = line + "#     Position in contig:\n"
    line = line + "#         title: 'Position in contig'\n"
    line = line + "#         description: 'Position in contig'\n"
    line = line + "#     Note:\n"
    line = line + "#         title: 'Note'\n"
    line = line + "#         description: 'Note'\n"
    line = line + "#     Accession number:\n"
    line = line + "#         title: 'Accession number'\n"
    line = line + "#         description: 'Accession number'\n"
    with open(plasmidfinder_output,'w') as f:
        f.write(line)

    plasmidfinder_df = pd.read_table(plasmidfinder_input)
    #plasmidfinder_df = plasmidfinder_df.iloc[:, :5]
    plasmidfinder_df.to_csv(plasmidfinder_output, sep="\t", mode='a')
    #if [ -f 'plasmidfinder_result.tsv' ]   ; then awk '{ print \$1 "_" NR "\t" \$2 "\t" \$3 "\t" \$4 "\t" \$5 }' plasmidfinder_result.tsv > plasmidfinder_mqc.tsv   ; fi

if exists(seqsero2_input):
    # sample  Output directory        Input files     O antigen prediction    H1 antigen prediction(fliC)     H2 antigen prediction(fljB)     Predicted identification        Predicted antigenic profile     Predicted serotype      Note
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'SeqSero2'\n"
    line = line + "# description: 'Salmonella serotype prediction from genome sequencing data.'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     O antigen prediction:\n"
    line = line + "#         title: 'O antigen prediction'\n"
    line = line + "#         description: 'O antigen prediction'\n"
    line = line + "#     H1 antigen prediction(fliC):\n"
    line = line + "#         title: 'H1 antigen prediction(fliC)'\n"
    line = line + "#         description: 'H1 antigen prediction(fliC)'\n"
    line = line + "#     H2 antigen prediction(fljB):\n"
    line = line + "#         title: 'H2 antigen prediction(fljB)'\n"
    line = line + "#         description: 'H2 antigen prediction(fljB)'\n"
    line = line + "#     Predicted identification:\n"
    line = line + "#         title: 'Predicted identification'\n"
    line = line + "#         description: 'Predicted identification'\n"
    line = line + "#     Predicted antigenic profile:\n"
    line = line + "#         title: 'Predicted antigenic profile'\n"
    line = line + "#         description: 'Predicted antigenic profile'\n"
    line = line + "#     Predicted serotype:\n"
    line = line + "#         title: 'Predicted serotype'\n"
    line = line + "#         description: 'Predicted serotype'\n"
    line = line + "#     Note:\n"
    line = line + "#         title: 'Note'\n"
    line = line + "#         description: 'Note'\n"
    with open(seqsero2_output,'w') as f:
        f.write(line)

    seqsero2_df = pd.read_table(seqsero2_input)
    #seqsero2_df = seqsero2_df.iloc[:, [0] + list(range(3, 10))]
    seqsero2_df.to_csv(seqsero2_output, index=False, sep="\t", mode='a')
    #if [ -f 'seqsero2_results.txt' ]       ; then cut -f 1,4-10 seqsero2_results.txt > seqsero2_mqc.txt ; fi

if exists(serotypefinder_input):
    # sample  Database        Gene    Serotype        Identity        Template / HSP length   Contig  Position in contig      Accession number
    line = "# plot_type: 'table'\n"
    line = line + "# section_name: 'SerotypeFinder'\n"
    line = line + "# description: 'SerotypeFinder identifies the serotype in total or partial sequenced isolates of E. coli.'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# headers:\n"
    line = line + "#     sample:\n"
    line = line + "#         title: 'sample'\n"
    line = line + "#         description: 'sample'\n"
    line = line + "#     Database:\n"
    line = line + "#         title: 'Database'\n"
    line = line + "#         description: 'Database'\n"
    line = line + "#     Gene:\n"
    line = line + "#         title: 'Gene'\n"
    line = line + "#         description: 'Gene'\n"
    line = line + "#     Serotype:\n"
    line = line + "#         title: 'Serotype'\n"
    line = line + "#         description: 'Serotype'\n"
    line = line + "#     Identity:\n"
    line = line + "#         title: 'Identity'\n"
    line = line + "#         description: 'Identity'\n"
    line = line + "#     Template / HSP length:\n"
    line = line + "#         title: 'Template / HSP length'\n"
    line = line + "#         description: 'Template / HSP length'\n"
    with open(serotypefinder_output,'w') as f:
        f.write(line)

    serotypefinder_df = pd.read_table(serotypefinder_input)
    serotypefinder_df = serotypefinder_df.iloc[:, :6]
    serotypefinder_df = serotypefinder_df.replace(to_replace=' ', value='', regex=True)
    serotypefinder_df.columns = serotypefinder_df.columns.str.replace(' ', '')
    serotypefinder_df.to_csv(serotypefinder_output, index=True, sep="\t", mode='a')
    #if [ -f 'serotypefinder_results.txt' ] ; then cut -f 1-6 serotypefinder_results.txt > serotypefinder_mqc.txt  ; fi
   
if exists(core_genome_input):
    # core,soft,shell,cloud
    line = "# plot_type: 'bargraph'\n"
    line = line + "# section_name: 'Core Genome Evaluation'\n"
    line = line + "# description: 'Genes Identified'\n"
    line = line + "# pconfig:\n"
    line = line + "#     namespace: 'Cust Data'\n"
    line = line + "# format: 'csv'\n"
    with open(core_genome_output,'w') as f:
        f.write(line)

    core_genome_df = pd.read_csv(core_genome_input)
    core_genome_df = core_genome_df.loc[:, ['sample', 'core', 'soft', 'shell', 'cloud']]
    core_genome_df.to_csv(core_genome_output, index=False, mode='a')


if exists(heatcluster_input):
    shutil.copyfile(heatcluster_input, heatcluster_output)

if exists(snpdists_input):
    shutil.copyfile(snpdists_input, snpdists_output)

if exists(phytreeviz_iqtr_input):
    shutil.copyfile(phytreeviz_iqtr_input, phytreeviz_iqtr_output)

if exists(phytreeviz_mshtr_input):
    shutil.copyfile(phytreeviz_mshtr_input, phytreeviz_mshtr_output)

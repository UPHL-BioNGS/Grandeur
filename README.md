# Grandeur

<img src="https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg" width="500" align="left" />

Named after [Grandeur Peak](https://www.alltrails.com/trail/us/utah/grandeur-peak-east-trail-from-church-fork).

**Image Credit:** [ryancornia](https://www.roadtripryan.com/go/resources/content/utah/wasatch/grandeur-peak/user-submitted/ryancornia-1505057017043.jpg)
**Location:** 40.707, -111.76, 8,299 ft (2,421 m) summit.
**Trail Info:** [https://utah.com/hiking/grandeur-peak](https://utah.com/hiking/grandeur-peak)

---

**Grandeur** is a species-agnostic sequencing analysis workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laboratory (UPHL)](https://uphl.utah.gov/). Built on [Nextflow](https://www.nextflow.io/), the pipeline provides quality control (QC), *de novo* assembly, taxonomic profiling, and in silico serotyping for paired-end Illumina data.

While intended to augment the CDC's [PHOENIX](https://github.com/CDCgov/phoenix) workflow, Grandeur also functions as a powerful standalone pipeline.

## Quick Start

### Dependencies
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (>= 25.0.0)
- [Apptainer/Singularity](https://apptainer.org/) or [Docker](https://www.docker.com/)

### Basic Usage
```bash
# Execution with FASTQ reads (Singularity)
nextflow run UPHL-BioNGS/Grandeur -profile singularity --reads <path_to_fastqs>

# Execution with existing Assemblies (Docker)
nextflow run UPHL-BioNGS/Grandeur -profile docker --fastas <path_to_fastas>

# Execution of a full Phylogenetic Analysis
nextflow run UPHL-BioNGS/Grandeur -profile singularity --sample_sheet samples.csv --msa
```
## Acknowledgements
## Workflow Architecture

Grandeur is modular and executes the following stages based on inputs and flags:
1.  **De Novo Alignment:** Cleaning reads with `fastp` and assembling with `SPAdes`.
2.  **Taxonomic Profiling:** Rapid identification via `SKANI`, `Kraken2`, `Mash`, and `Sylph`.
3.  **Quality Assessment:** Assembly metrics via `QUAST`, `CheckM2`, and `FastQC`.
4.  **Subtyping:** Organism-specific serotyping (e.g., `Kleborate`, `SeqSero2`, `Legionella` SBT).
5.  **Phylogenetic Analysis (Optional):** Core genome alignment and Maximum Likelihood trees.

---

## Key Parameters

| Parameter | Description |
| :--- | :--- |
| `--sample_sheet` | CSV with `sample,fastq_1,fastq_2` |
| `--reads` | Directory with paired-end FASTQ files |
| `--fastas` | Directory with FASTA files |
| `--outdir` | Directory to save results (default: `grandeur`) |
| `--msa` | **Toggle:** Run phylogenetic analysis |
| `--skip_extras` | **Toggle:** Run only core assembly/QC |
| `--current_datasets`| **Toggle:** Download NCBI references via `datasets` |

*The full list of parameters is available by running:* `nextflow run UPHL-BioNGS/Grandeur --help`

---

## Documentation
Detailed guides, FAQ, and process explanations are located in the **[Grandeur Wiki](https://github.com/UPHL-BioNGS/Grandeur/wiki)**.

- [Installation](https://github.com/UPHL-BioNGS/Grandeur/wiki/Installation)
- [Usage Examples](https://github.com/UPHL-BioNGS/Grandeur/wiki/examples)
- [Phylogenetic Analysis](https://github.com/UPHL-BioNGS/Grandeur/wiki/Phylogenetic-Analysis)
- [Subworkflow Explanations](https://github.com/UPHL-BioNGS/Grandeur/wiki/subworkflows)

## Acknowledgements
Grandeur integrates an extensive suite of open-source bioinformatics tools, including:
`amrfinderplus`, `bakta`, `checkm2`, `drprg`, `elgato`, `emmtyper`, `fastp`, `fastqc`, `gotree`, `heatcluster`, `iqtree`, `kaptive`, `kleborate`, `kraken2`, `ksnp4`, `mash`, `mashtree`, `meningotype`, `mlst`, `multiqc`, `mykrobe`, `ngmaster`, `panaroo`, `pbptyper`, `plasmidfinder`, `prokka`, `quast`, `roary`, `seqsero2`, `serotypefinder`, `shigapass`, `ska2`, `skani`, `snp-dists`, `spades`, `sylph`.

## Technical Support
Issues and problems should be submitted to the [GitHub Issues](https://github.com/UPHL-BioNGS/Grandeur/issues) page.

Grandeur wouldn't be possible without the following tools:
- nf-tools
- [amrfinderplus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) - identification of genes associated with antimicrobial resistence
- bakta
- checkm2
- [datasets](https://github.com/ncbi/datasets) - downloads genomes from NCBI
- [drprg](https://github.com/mbhall88/drprg) - TB AMR predictions
- [elgato](https://github.com/appliedbinf/el_gato) - Legionella pneumophila Sequence Based Typing (SBT)
- [emmtyper](https://github.com/MDU-PHL/emmtyper) - Group A Strep "emm" typing
- enatools
- [fastp](https://github.com/OpenGene/fastp) - cleaning reads
- [fastqc](https://github.com/s-andrews/FastQC) - fastq file QC
- gotree
- [heatcluster](https://github.com/erinyoung/heatcluster) - visualizes SNP matrix from SNP dists
- [iqtree](http://www.iqtree.org/) - phylogenetic tree creation - used after core genome alignment
- kaptive
- [kleborate](https://github.com/katholt/Kleborate) - Klebsiella serotyping
- [kraken2](https://ccb.jhu.edu/software/kraken2/) - contamination
- ksnp4
- [mash](https://github.com/marbl/Mash) - species identifier
- [mashtree](https://github.com/lskatz/mashtree) - tree based on mash distances (not impacted by size of core genome)
- meningotype
- [mlst](https://github.com/tseemann/mlst) - identification of MLST subtype
- [multiqc](https://multiqc.info/) - summarizes QC efforts
- [mykrobe](https://github.com/Mykrobe-tools/mykrobe) - Mycobacterium subtyping
- ngmaster
- [panaroo](https://github.com/gtonkinhill/panaroo) - core genome alignment - optional (set with params.msa = true)
- [pbptyper](https://github.com/rpetit3/pbptyper) - Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae assemblies
- [plasmidfinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/) - MLST typing for plasmids
- [prokka](https://github.com/tseemann/prokka) - gene annotation - used for core genome alignment
- [quast](http://quast.sourceforge.net/quast) - contig QC
- roary
- [seqsero2](https://github.com/denglab/SeqSero2) - Salmonella serotyping
- [seqsero2S]() - Salmonella serotyping
- [serotypefinder](https://cge.cbs.dtu.dk/services/SerotypeFinder/) - E. coli serotyping
- [shigapass]() - Shigella serotyping
- ska2
- skani
- [snp-dists](https://github.com/tseemann/snp-dists) - SNP matrix - used after core genome aligment
- [spades](https://cab.spbu.ru/software/spades/) - _de novo_ alignment
- spestimator
- sylph

The expected tools are split into multiple processes. Each [process has its own wiki page](https://github.com/UPHL-BioNGS/Grandeur/wiki/Processes) that we encourage users to view.

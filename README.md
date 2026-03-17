# Sponge genome decontamination workflow

This repository contains the scripts used in the analysis presented in the paper: “Hidden Contaminants in Sponge Genomes: Large-Scale Decontamination of 30 Public Assemblies” (Bodulić and Vlahoviček, 2026). More information can be found in the preprint (coming soon)

# General Information

To run the full pipeline, clone this repository and execute the master script (sponge_decontamination_workflow.bash) within the repository environment:

```bash
./sponge_decontamination_workflow.bash
```
Make sure all required dependencies are added to the `PATH` environment variable. All scripts in the scripts directory should also be added to `PATH` before running the script. An active internet connection is required.

# Dependencies

The following dependencies were used in the pipeline:

| **Dependency**                   | **Version** |                                                              
|----------------------------------|-------------|
| datasets                         | 18.20.0     |
| R                                | 4.5.2       |
| stringr (R)                      | 1.5.1       | 
| scales  (R)                      | 1.4.0       |
| Genome decontamination pipeline  | 1.0.0       |
| BUSCO                            | 6.0.0       |
| wesanderson (R)                  | 0.3.7       |
| ComplexUpset (R)                 | 1.3.3       |
| pheatmap (R)                     | 1.0.13      |
| taxsize (R)                      | 0.9.98      |
| ggVennDiagram (R)                | 1.5.4       |
| patchwork (R)                    | 1.3.2       |

Tools denoted with (R) correspond to R packages. Version information for dependencies specific to genome decontamination pipeline can be found in the README file of its repository.

The workflow was run on a high-performance computing cluster using SUSE Linux Enterprise High Performance Computing 15 SP4 (SLE_HPC 15-SP4) with Linux kernel 5.14.21. The system used a GNU Bash shell version 4.4.23, and ran on a node equipped with Intel(R) Xeon(R) Gold 6248 CPUs @ 2.50GHz. 

# Script overview

This repository contains the following scripts:

| **Script**                             | **Purpose**                                                 |        
|----------------------------------------|-------------------------------------------------------------|
| `sponge_decontamination_workflow.bash` | Master script for the analysis                              |
| `generate_validation_dataset.R`        | Generates validation spike-in dataset from source sequences |
| `generate_fig2_elements.R`             | Generates Figure 2 elements                                 |
| `generate_fig3_elements.R`             | Generates Figure 3 elements                                 |

The repository also contains the following files:

| **File**                       | **Purpose**                                                                                   |        
|--------------------------------|-----------------------------------------------------------------------------------------------|
| `analysed_sponge_genomes.tsv`  | Contains accession numbers of the analysed sponge genome assemblies (NCBI Genomes)            |
| `taxonomy_info_validation.csv` | Contains full taxonomic classification of the contaminants included in the validation dataset |

For convenience, tables used in figure generation are also directly provided:

| **Table**                                         | **Figure**          | **Download link**                                                                                                 | 
|---------------------------------------------------|---------------------|-------------------------------------------------------------------------------------------------------------------|
| `scaffold_contamination_info_validation.tsv`      | Figure 2 (a-c)      | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/scaffold_contamination_info_validation.tsv.gz)      |
| `contaminant_scaffolds_taxonomy_validation.tsv`   | Figure 2 (d)        | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/contaminant_scaffolds_taxonomy_validation.tsv.gz)   |
| `taxonomy_info_validation.csv`                    | Figure 2 (c,d)      | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/taxonomy_info_validation.csv.gz)                    |
| `scaffold_contamination_info_sponges.tsv`         | Figure 3 (a-c, e,f) | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/scaffold_contamination_info_sponges.tsv.gz)         |
| `contaminant_scaffolds_taxonomy_sponges.tsv`      | Figure 3 (d)        | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/contaminant_scaffolds_taxonomy_sponges.tsv.gz)      |
| `mc_umap_dt_sponges.tsv`                          | Figure 3 (f)        | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/mc_umap_dt_sponges.tsv.gz)                          |
| `busco_original_genome_complete_merged.tsv`       | Figure 3 (g)        | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/busco_original_genome_complete_merged.tsv.gz)       |
| `busco_decontaminated_genome_complete_merged.tsv` | Figure 3 (g)        | [Download](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/busco_decontaminated_genome_complete_merged.tsv.gz) |

Supplementary Figure 1 is automatically generated when creating the validation spike-in dataset (script `generate_validation_dataset.R`)

Tables should be unzipped before being directly supplied to the corresponding R scripts.

If the files are not directly accessible via the provided links, they are available in the following [folder](http://hex.bioinfo.hr/~kbodulic/genom_decom_tables/).

The tables have also been deposited to Zenodo (10.5281/zenodo.19064896)

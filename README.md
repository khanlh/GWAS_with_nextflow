# GWAS Nextflow Demo

This repository contains a **minimal reproducible workflow** for performing a **Genome-Wide Association Study (GWAS)** using [Nextflow](https://www.nextflow.io/) and Docker.  
It uses a **mini simulated dataset** so that the workflow can be tested quickly without requiring large genomic data.

---

## 1. What is GWAS?

Genome-Wide Association Studies (GWAS) are a method to identify associations between genetic variants (usually SNPs) and traits (phenotypes).  
The general steps are:

1. **Genotype Data Preparation**  
   - SNPs are genotyped for many individuals, usually stored in PLINK binary format (`.bed`, `.bim`, `.fam`).  

2. **Phenotype Collection**  
   - Phenotype data (e.g., disease status, height, weight) is prepared in a simple text file.  

3. **Quality Control (QC)**  
   - Remove low-quality SNPs or individuals.  
   - Filter by call rate, Hardy-Weinberg equilibrium, or minor allele frequency.  

4. **Association Testing**  
   - Statistical tests are performed for each SNP against the phenotype.  
   - Tools like **PLINK** calculate p-values for SNP-trait associations.  

5. **Visualization**  
   - Results are visualized in a **Manhattan plot** (−log10(p-value) vs genomic position).  
   - A **QQ plot** can be used to evaluate the distribution of p-values.  

---

## 2. Workflow Overview

This demo workflow implements a simplified GWAS pipeline:

- **Step 1: Simulate or provide genotype + phenotype data**  
  - We use `plink --simulate-geno` to create a small test dataset.  

- **Step 2: Association Testing**  
  - Run `plink --assoc` with the genotype and phenotype files.  

- **Step 3: Visualization**  
  - Plot a Manhattan plot using R.  

---

## 3. Repository Structure
gwas-nextflow-demo/
├── data/
│ ├── genotype.bed # Genotype binary file
│ ├── genotype.bim # SNP information
│ ├── genotype.fam # Individual information
│ └── phenotype.txt # Phenotype data
├── Dockerfile # Environment with PLINK + R
├── main.nf # Nextflow workflow
├── nextflow.config # Config (Docker support)
└── README.md # This file


---

## 4. Requirements

- [Nextflow](https://www.nextflow.io/) (≥22.10)  
- [Docker](https://docs.docker.com/)  

---

## 5. Build Docker Image

```bash
docker build -t gwas-demo:latest .


## 6. Run the Workflow
nextflow run main.nf -with-docker gwas-demo:latest

## 7. Output Files

gwas_results.assoc → Raw GWAS association results from PLINK

manhattan.png → A simple Manhattan plot generated from the results

## 8. Example Manhattan Plot

(For illustration, real output will depend on the simulated dataset)

## 9. Notes

This is a minimal demo, not a production-ready GWAS pipeline.

For real GWAS projects:

Perform extensive QC (SNP/individual filters).

Use covariates (age, sex, population structure).

Consider mixed models (e.g., GEMMA, GCTA).

## 10. Workflow Diagram

```mermaid
flowchart TD
    A[Genotype + Phenotype Data] --> B[Association Testing (PLINK)]
    B --> C[Manhattan Plot Visualization (R)]


---
##11. References

Purcell S, et al. (2007) PLINK: A tool set for whole-genome association and population-based linkage analyses.

Tutorial on GWAS: https://www.ebi.ac.uk/training/online/courses/gwas-tutorial/

Nextflow documentation: https://www.nextflow.io/docs/latest/index.html

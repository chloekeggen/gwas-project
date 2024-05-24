# gwas-project
# Genome-Wide Association Study (GWAS) Pipeline

This project implements a pipeline to perform a Genome-Wide Association Study (GWAS) using Python. The pipeline includes steps for data loading, quality control, Hardy-Weinberg Equilibrium (HWE) testing, filtering, statistical testing, and result visualization.

## Project Overview

The pipeline performs the following tasks:
1. Data Import and Management
2. Quality Control
3. Association Testing
4. Result Visualization
5. Multiple Testing Correction

## Setup Instructions

### Prerequisites

Ensure you have the following software installed:
- Python 3.6+
- Pip package manager

### Install Required Python Packages

Install the required Python packages using pip: ```bash
pip install pandas numpy matplotlib cyvcf2 scipy


## Data Files
Place the following data files in the project directory:

lab3_gwas.phen: Phenotype file with normalized LDL values.
lab3_gwas.vcf.gz: VCF f##ile containing LD-pruned SNPs.
lab3_gwas.vcf.gz.tbi: Index file for the VCF file.

##Data Import and Management
### Functions
read_geno_vcf(file_path): Reads VCF file and returns genotype data as a pandas DataFrame.
read_pheno_file(file_path): Reads phenotype data from a CSV file and returns it as a pandas DataFrame.

## Quality Control
### Functions
filter_snps(genotype_df, maf_threshold, missingness_threshold): Filters SNPs based on minor allele frequency (MAF) and missingness thresholds.
filter_indivs(genotype_df, missingness_threshold): Filters individuals based on a missingness threshold.
filter_hwe(genotype_df, hwe_threshold=1e-6): Filters SNPs based on Hardy-Weinberg Equilibrium (HWE) p-values.

## Association Testing
### Functions
linear_regression(genotype_df, phenotype_df): Performs linear regression to find associations between SNPs and phenotypes.

## Result Visualization
### Functions
manhattan_plot(results): Generates a Manhattan plot to visualize GWAS results.
qq_plot(results): Generates a Q-Q plot to depict the distribution of p-values.


###Running the Analysis
Load Phenotype Data: Use read_pheno_file to load phenotype data.
Load Genotype Data: Use read_geno_vcf to load genotype data from the VCF file.
Quality Control:
Use filter_snps to filter SNPs based on MAF and missingness.
Use filter_indivs to filter individuals based on missingness.
Use filter_hwe to filter SNPs based on HWE.
Association Testing: Use linear_regression to find associations between SNPs and phenotypes.
Visualization: Use manhattan_plot and qq_plot to visualize the results.

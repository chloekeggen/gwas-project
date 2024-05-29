# Linear GWAS Package

This package provides a simple tool for performing Genome-Wide Association Studies (GWAS) on continuous phenotypes using linear regression.

## Installation

To install the package, access via GitHub link:

```
git clone https://github.com/chloekeggen/gwas-project.git
cd gwas-project
python setup.py install
```

## Usage

After installing the package, you can use gwas-tools-cli.py to perform GWAS on your data.

```
gwas-tools-cli --vcf <path_to_vcf_file> --pheno <path_to_phenotype_file> --out <output_file_prefix>
```

Replace <path_to_vcf_file> with the path to your VCF file containing genotype data, <path_to_phenotype_file> with the path to your phenotype file, and <output_file_prefix> with the desired name for the output files.

### Optional arguments

```
--maf <maf_threshold>
```

Adjust MAF threshold for filtering SNPs as needed; the default is 0.05.

## Output

- <output_file_prefix>_results.csv: CSV file containing the results of the linear regression analysis, including SNP IDs, chromosome, position, and p-values.
- <output_file_prefix>_manhattan_plot.png: Manhattan plot visualizing the results of the GWAS analysis.

## Example using given phenotype and genotype files

```
gwas-tools-cli --vcf https://github.com/chloekeggen/gwas-project/raw/main/lab3_gwas.vcf.gz --pheno https://github.com/chloekeggen/gwas-project/raw/main/lab3_gwas.phen --out gwas_results
```

This command will perform GWAS on the provided genotype and phenotype data files from Lab 3, and save the results into 2 gwas_results files.

## Dependencies
- pandas
- numpy
- cyvcf2
- scikit-learn
- statsmodels
- matplotlib

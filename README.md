# GWAS Package

A package for preprocessing GWAS data.

## Installation

```sh
pip install git+https://github.com/chloekeggen/gwas-project.git

USAGE
from gwas.preprocessing import preprocess_gwas_data

# Replace these file paths with the actual paths to your phenotype and genotype files
phenotype_file = "path/to/phenotype_file"
genotype_file = "path/to/genotype_file"

# Preprocess the GWAS data
phenotype_data, genotype_data = preprocess_gwas_data(phenotype_file, genotype_file)

# Print the shapes of the processed data
print("Genotype data info:", genotype_data.shape)
print("Phenotype data info:", phenotype_data.shape)


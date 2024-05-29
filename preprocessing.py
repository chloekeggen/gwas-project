import os
import pandas as pd
import numpy as np
from cyvcf2 import VCF

def read_phenotype_data(file_path):
    return pd.read_csv(file_path, delim_whitespace=True, names=["FID", "IID", "LDL"])

def read_genotype_data(file_path):
    samples = (VCF(file_path)).samples
    genotype_dict = {}

    for data in VCF(file_path):
        genotype_data = [gt[0] + gt[1] for gt in data.genotypes]
        genotype_dict[data.ID] = genotype_data

    genotype_data = pd.DataFrame(genotype_dict, index=samples)
    return genotype_data

def filter_snps(genotype_df, maf_threshold=0.05, missingness_proportion=0.05):
    maf = genotype_df.mean(axis=0) / 2
    snps = maf[maf > maf_threshold].index
    missingness_prop = genotype_df.isnull().mean()
    snps = snps.intersection(missingness_prop[missingness_prop < missingness_proportion].index)
    return genotype_df[snps]

def preprocess_gwas_data(phenotype_file, genotype_file):
    phenotype_data = read_phenotype_data(phenotype_file)
    genotype_data = read_genotype_data(genotype_file)

    # Perform QC
    genotype_data = filter_snps(genotype_data)

    return phenotype_data, genotype_data

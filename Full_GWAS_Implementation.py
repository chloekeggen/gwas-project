#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
from cyvcf2 import VCF
from sklearn.decomposition import PCA
import statsmodels.api as sm
import matplotlib.pyplot as plt

def read_phenotype_data(file_path):
    return pd.read_csv(file_path, delim_whitespace=True, names=["FID", "IID", "LDL"])

def read_genotype_data(file_path, chromosomes=None):
    samples = (VCF(file_path)).samples
    genotype_dict = {}
    chrom_dict = {}
    pos_dict = {}

    for data in VCF(file_path):
        if chromosomes is None or data.CHROM in chromosomes:
            genotype_data = [gt[0] + gt[1] for gt in data.genotypes]
            genotype_dict[data.ID] = genotype_data
            chrom_dict[data.ID] = data.CHROM
            pos_dict[data.ID] = data.POS

    genotype_data = pd.DataFrame(genotype_dict, index=samples)
    return genotype_data, chrom_dict, pos_dict


# In[2]:


# Filter SNPs based on criteria used by the Wellcome Trust Case Control Consortium (WTCCC) when conducting GWAS
# The criteria for retaining a SNP are: HWE P-value ≥ 5.7 × 10−7, MSP ≤5% if MAF ≥ 5%, MSP ≤ 1%
def filter_snps(genotype_df, maf_threshold = 0.05, missingness_proportion = 0.05):
    maf = genotype_df.mean(axis = 0) / 2
    snps = maf[maf > maf_threshold].index
    missingness_prop = genotype_df.isnull().mean()
    snps = snps.intersection(missingness_prop[missingness_prop < missingness_proportion].index)
    return genotype_df[snps]

# Follow same guidelines usde by WTCCC to filter individuals based on missingness
def filter_indivs(genotype_df, missingness_proportion = 0.05):
    missingness = genotype_df.isnull().mean(axis = 1)
    genotype_df = genotype_df.loc[missingness < missingness_proportion, :]
    return genotype_df


# In[3]:


# Main function to perform GWAS data preprocessing
def preprocess_gwas_data(phenotype_file, genotype_file):
    phenotype_data = read_phenotype_data(phenotype_file)
    chromosomes = [str(i) for i in range(1, 3)]
    genotype_data, chrom_dict, pos_dict = read_genotype_data(genotype_file, chromosomes=chromosomes)
    
    # Perform QC
    genotype_data = filter_snps(genotype_data)
    genotype_data = filter_indivs(genotype_data)
    
    return phenotype_data, genotype_data, chrom_dict, pos_dict


# In[4]:


def perform_pca(genotype_data, n_components=3):
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(genotype_data.fillna(0))
    pcs_df = pd.DataFrame(pcs, index=genotype_data.index, columns=[f'PC{i+1}' for i in range(n_components)])
    return pcs_df


# In[5]:


def perform_linear_regression(phenotype_data, genotype_data, pcs_df, chrom_dict, pos_dict):
    phenotype_data.reset_index(drop=True, inplace=True)
    genotype_data.reset_index(drop=True, inplace=True)
    pcs_df.reset_index(drop=True, inplace=True)

    phenotype_data = phenotype_data.apply(pd.to_numeric, errors='coerce')
    y = phenotype_data['LDL']

    results = []

    for snp_col in genotype_data.columns:
        if snp_col in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
            continue
            
        chrom = chrom_dict[snp_col]
        pos = pos_dict[snp_col]

        genotype_info = genotype_data[snp_col]

        if genotype_info.isnull().any():
            print(f"Genotype data for SNP {snp_col} contains missing values. Skipping this SNP.")
            continue

        X = pcs_df.copy()
        X.insert(0, 'const', 1)
        X[snp_col] = genotype_info

        # Peforming linear regression
        model = sm.OLS(y, X).fit()
        p_value = model.pvalues[snp_col]
        results.append((snp_col, chrom, pos, p_value))

    return results


# In[6]:


def generate_manhattan_plot(phenotype_file, genotype_file, output_file):
    phenotype_file = os.path.expanduser(phenotype_file)
    genotype_file = os.path.expanduser(genotype_file)
    output_file = os.path.expanduser(output_file)
    
    phenotype_data, genotype_data, chrom_dict, pos_dict = preprocess_gwas_data(phenotype_file, genotype_file)
    pcs_df = perform_pca(genotype_data, n_components=3)
    results = perform_linear_regression(phenotype_data, genotype_data, pcs_df, chrom_dict, pos_dict)

    results_df = pd.DataFrame(results, columns=['SNP', 'Chromosome', 'Position', 'P-Value'])
    results_df.to_csv("gwas_results.csv", index=False)

    # Calculating position on chromosome for plotting
    results_df['Chromosome'] = results_df['Chromosome'].astype(int)
    results_df = results_df.sort_values(['Chromosome', 'Position'])
    chrom_max = results_df.groupby('Chromosome')['Position'].max().cumsum()
    results_df['Cumulative_Position'] = results_df.apply(lambda row: row['Position'] + (chrom_max[row['Chromosome']] if row['Chromosome'] in chrom_max else 0), axis=1)

    # Creating Manhattan plot
    plt.figure(figsize=(15, 5))
    for chrom, group in results_df.groupby('Chromosome'):
        color = 'black' if chrom % 2 == 0 else 'gray'
        plt.scatter(group['Cumulative_Position'], -np.log10(group['P-Value']), c=color, s=10)

    plt.axhline(-np.log10(0.0001), color='blue', linestyle='--', linewidth=1)
    plt.axhline(-np.log10(0.00001), color='red', linestyle='--', linewidth=1)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(P-Value)')
    plt.title('Manhattan Plot')

    chrom_ticks = results_df.groupby('Chromosome')['Cumulative_Position'].median()
    plt.xticks(chrom_ticks, chrom_ticks.index)

    plt.savefig(output_file)
    plt.close()


# In[7]:


# Accepting user input for file paths
phenotype_file = input("Enter the path to the phenotype file: ")
genotype_file = input("Enter the path to the genotype file: ")
output_file = input("Enter the path for saving the Manhattan plot: ")


# In[8]:


generate_manhattan_plot(phenotype_file, genotype_file, output_file)


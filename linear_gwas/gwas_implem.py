import os
import pandas as pd
import numpy as np
import vcf
import statsmodels.api as sm
import matplotlib.pyplot as plt
import warnings

# Suppressing FutureWarning from statsmodels
warnings.filterwarnings("ignore", category = FutureWarning)

def read_phenotype_data(file_path):
    """
    Read phenotype data from CSV file.

    :param file_path: Path to the phenotype file.
    :return: Phenotype data as pandas.DataFrame.
    """
    return pd.read_csv(file_path, delim_whitespace = True, names = ["FID", "IID", "LDL"])

def read_genotype_data(file_path):
    """
    Read genotype data from a VCF file.

    :param file_path: Path to the VCF file.
    :return: Genotype data as pandas.DataFrame,
             Dictionary with chromosome info for each SNP,
             Dictionary containing position on the chromosome info for each SNP.
    """
    vcf_reader = vcf.Reader(filename = file_path)
    samples = vcf_reader.samples
    genotype_dict = {}
    chrom_dict = {}
    pos_dict = {}

    for record in vcf_reader:
        genotype_data = [sum(map(int, call.data.GT.split('|'))) for call in record.samples]
        genotype_dict[record.ID] = genotype_data
        chrom_dict[record.ID] = record.CHROM
        pos_dict[record.ID] = record.POS

    genotype_data = pd.DataFrame(genotype_dict, index = samples)
    
    return genotype_data, chrom_dict, pos_dict


def compute_maf_plink_like(genotype_SNP):
    """
    Compute Minor Allele Frequency (MAF) for a genotype series.

    :param genotypes: Genotype data for SNP.
    :return: MAF.
    """
    # Counting number of occurrences of each genotype (0, 1, 2), assuming SNPs are bi-allelic
    genotype_counts = genotype_SNP.value_counts(dropna=True)
    allele_counts = {0: 0, 1: 0}
    total_alleles = 0

    # Counting alleles
    for genotype, count in genotype_counts.items():
        if genotype == 0:
            allele_counts[0] += 2 * count
        elif genotype == 1:
            allele_counts[0] += count
            allele_counts[1] += count
        elif genotype == 2:
            allele_counts[1] += 2 * count
        total_alleles += 2 * count

    # Calculating MAF
    if total_alleles == 0:
        return 0.0
    
    minor_allele_count = min(allele_counts.values())
    maf = minor_allele_count / total_alleles

    return maf

def filter_snps(genotype_df, maf_threshold=0.05):
    """
    Filtering SNPs based on MAF threshold.

    :param genotype_df: Genotype dataframe.
    :param maf_threshold: Minimum MAF threshold.
    :return: Filtered genotype data based on MAF threshold.
    """
    mafs = genotype_df.apply(compute_maf_plink_like)
    snps = mafs[mafs >= maf_threshold].index
    return genotype_df[snps]

def preprocess_gwas_data(phenotype_file, genotype_file):
    """
    Preprocessing GWAS data.

    :param phenotype_file: Path to the phenotype file.
    :param genotype_file: Path to the genotype file.
    :return: Phenotype data, genotype data, chromosome information dictionary, position information dictionary.
    """
    phenotype_data = read_phenotype_data(phenotype_file)
    genotype_data, chrom_dict, pos_dict = read_genotype_data(genotype_file)
    
    genotype_data = filter_snps(genotype_data)
    
    return phenotype_data, genotype_data, chrom_dict, pos_dict

def merge_genotype_phenotype(phenotype_data, genotype_data):
    """
    Merging phenotype and genotype data.

    :param phenotype_data: Phenotype data.
    :param genotype_data: Genotype data.
    :return: Merged phenotype and genotype data.
    """
    merged_data = pd.merge(phenotype_data, genotype_data, left_on='IID', right_index=True)
    return merged_data

def perform_linear_regression(phenotype_data, genotype_data, chrom_dict, pos_dict):
    """
    Performing linear regression for SNP.

    :param phenotype_data: Phenotype data.
    :param genotype_data: Genotype data.
    :param chrom_dict: Chromosome information dictionary.
    :param pos_dict: Position information dictionary.
    :return: List of linear regression results.
    """
    merged_data = merge_genotype_phenotype(phenotype_data, genotype_data)
    
    y = merged_data['LDL']
    results = []

    for snp_col in genotype_data.columns:
        if snp_col in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
            continue
            
        chrom = chrom_dict[snp_col]
        pos = pos_dict[snp_col]

        genotype_info = merged_data[snp_col]

        if genotype_info.isnull().any():
            print(f"Genotype data for SNP {snp_col} has missing info, so SNP skipped.")
            continue

        # Assuming A2 (major allele) by PLINK convention, therefore flip genotypes, so A1 (minor allele) is the reference
        genotype_info = 2 - genotype_info
        
        # Add intercept manually
        X = sm.add_constant(genotype_info)
        
        model = sm.OLS(y, X).fit()
        p_value = model.pvalues[1]
        beta = model.params[1]
        t_stat = model.tvalues[1]
        
        results.append((chrom, snp_col, pos, 'ADD', len(y), beta, t_stat, p_value))

    return results
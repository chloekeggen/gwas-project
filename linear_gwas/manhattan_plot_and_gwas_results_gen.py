import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from .gwas_implem import preprocess_gwas_data, perform_linear_regression

def generate_manhattan_plot_and_gwas_results(results_df, output_prefix):
    """
    Generate a Manhattan plot for GWAS.

    :param results_df: DataFrame containing GWAS results.
    :param output_file: Path to save Manhattan plot.
    """
    output_prefix = os.path.expanduser(output_prefix)
    
    results_df['CHR'] = results_df['CHR'].astype(int)
    results_df = results_df.sort_values(['CHR', 'BP'])
    chrom_max = results_df.groupby('CHR')['BP'].max().cumsum()
    results_df['final_post'] = results_df.apply(lambda row: row['BP'] + (chrom_max[row['CHR']] if row['CHR'] in chrom_max else 0), axis = 1)
    
    plt.figure(figsize=(15, 5))
    for chrom, group in results_df.groupby('CHR'):
        color = 'black' if chrom % 2 == 0 else 'gray'
        plt.scatter(group['final_post'], -np.log10(group['P']), c = color, s = 10)

    plt.axhline(-np.log10(0.00001), color = 'blue', linestyle = '--', linewidth = 1)
    plt.axhline(-np.log10(0.0000001), color = 'red', linestyle = '--', linewidth = 1)
    plt.xlabel('Chromosomes')
    plt.ylabel('-log10(p)')
    plt.title('Manhattan Plot')

    chrom_ticks = results_df.groupby('CHR')['final_post'].median()
    plt.xticks(chrom_ticks, chrom_ticks.index)

    plt.savefig(f"{output_prefix}_manhattan_plot.png")
    plt.close()
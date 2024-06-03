import numpy as np
import matplotlib.pyplot as plt
import os

def generate_qq_plot(p_values, output_prefix):
    """
    Generate a QQ plot for GWAS.

    :param p_values: Array of p-values.
    :param output_prefix: Prefix for the output files.
    """
    output_prefix = os.path.expanduser(output_prefix)
    
    # Sort the p-values
    sorted_p_values = np.sort(p_values)
    expected_p_values = np.arange(1, len(sorted_p_values) + 1) / len(sorted_p_values)
    
    plt.figure(figsize=(6, 6))
    plt.scatter(-np.log10(expected_p_values), -np.log10(sorted_p_values), c='blue', s=10)
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlabel('Expected -log10(p)')
    plt.ylabel('Observed -log10(p)')
    plt.title('QQ Plot')
    
    plt.savefig(f"{output_prefix}_qq_plot.png")
    plt.close()
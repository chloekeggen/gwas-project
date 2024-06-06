import argparse
import pandas as pd
import os
from .gwas_implem import preprocess_gwas_data, perform_linear_regression
from .manhattan_plot_and_gwas_results_gen import generate_manhattan_plot_and_gwas_results
from .qq_plot import generate_qq_plot

def main():
    parser = argparse.ArgumentParser(description = "Simple GWAS on continuous phenotypes via linear regression")
    
    parser.add_argument("--vcf", required = True, help = "Path to genoptype VCF file")
    parser.add_argument("--maf", type = float, default = 0.05, help = "MAF threshold for filtering SNPs")
    parser.add_argument("--pheno", required = True, help = "Path to phenotype CSV file")
    parser.add_argument("--out", required = True, help = "Output file prefix")
    
    args = parser.parse_args()

    # Check if user input is correct
    for file_path in [args.vcf, args.pheno]:
        if not os.path.isfile(file_path):
            print(f"Please provide a correct file path. File not found: {file_path}")
            return
    
    phenotype_file = args.pheno
    genotype_file = args.vcf
    maf_threshold = args.maf
    output_prefix = args.out

    try:
        phenotype_data, genotype_data, chrom_dict, pos_dict = preprocess_gwas_data(phenotype_file, genotype_file)
    except Exception as e:
        print(f"Issue when preprocessing GWAS data: {e}")
        return
        
    try:
        results = perform_linear_regression(phenotype_data, genotype_data, chrom_dict, pos_dict)
    except Exception as e:
        print(f"Issue when performing linear regression: {e}")
        return

    results_df = pd.DataFrame(results, columns = ['CHR', 'SNP', 'BP', 'TEST', 'NMISS', 'BETA', 'STAT', "P"])
    results_df.to_csv(f"{output_prefix}_results.csv", index = False)
    
    generate_manhattan_plot_and_gwas_results(results_df, f"{output_prefix}_manhattan_plot.png")
    generate_qq_plot(results_df['P'], output_prefix)

    print(f"Results successfully saved to {output_prefix}_gwas_results.csv")
    print(f"Manhattan plot successfully saved to {output_prefix}_manhattan_plot.png")
    print(f"QQ plot successfully saved to {output_prefix}_qq_plot.png")

if __name__ == "__main__":
    main()

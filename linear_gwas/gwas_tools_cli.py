import argparse
import pandas as pd
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
    
    phenotype_file = args.pheno
    genotype_file = args.vcf
    maf_threshold = args.maf
    output_prefix = args.out

    phenotype_data, genotype_data, chrom_dict, pos_dict = preprocess_gwas_data(phenotype_file, genotype_file)
    
    results = perform_linear_regression(phenotype_data, genotype_data, chrom_dict, pos_dict)
    results_df = pd.DataFrame(results, columns = ['CHR', 'SNP', 'BP', 'TEST', 'NMISS', 'BETA', 'STAT', "P"])
    results_df.to_csv(f"{output_prefix}_results.csv", index = False)
    
    generate_manhattan_plot_and_gwas_results(results_df, f"{output_prefix}_manhattan_plot.png")
    generate_qq_plot(results_df['P'], output_prefix)

    print(f"Results successfully saved to {output_prefix}_gwas_results.csv")
    print(f"Manhattan plot successfully saved to {output_prefix}_manhattan_plot.png")
    print(f"QQ plot successfully saved to {output_prefix}_qq_plot.png")

if __name__ == "__main__":
    main()
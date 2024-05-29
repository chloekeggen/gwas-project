import argparse
import pandas as pd
from .other_modules import preprocess_gwas_data, perform_pca, perform_linear_regression

def main():
    parser = argparse.ArgumentParser(description="Simple GWAS on continuous phenotypes via linear regression")
    
    parser.add_argument("--vcf", required=True, help="Path to VCF file")
    parser.add_argument("--maf", type=float, default=0.05, help="MAF threshold for filtering SNPs")
    parser.add_argument("--pheno", required=True, help="Path to phenotype file")
    parser.add_argument("--out", required=True, help="Output file prefix")
    
    args = parser.parse_args()
    
    phenotype_file = args.pheno
    genotype_file = args.vcf
    maf_threshold = args.maf
    output_prefix = args.out

    phenotype_data, genotype_data, chrom_dict, pos_dict = preprocess_gwas_data(phenotype_file, genotype_file, maf_threshold)
    
    pcs_df = perform_pca(genotype_data, n_components=3)
    results = perform_linear_regression(phenotype_data, genotype_data, pcs_df, chrom_dict, pos_dict)
    
    results_df = pd.DataFrame(results, columns=['SNP', 'Chromosome', 'Position', 'P-Value'])
    results_df.to_csv(f"{output_prefix}_results.csv", index=False)
    
    generate_manhattan_plot(results_df, output_prefix)

    print(f"Results saved to {output_prefix}_results.csv")
    print(f"Manhattan plot saved to {output_prefix}_manhattan_plot.png")

if __name__ == "__main__":
    main()
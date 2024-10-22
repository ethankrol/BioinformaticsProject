import pandas as pd

num_genes = [10, 100, 1000, 10000]

for num_genes in num_genes:
    # load the filtered data
    df = pd.read_csv("results/filtered_data.tsv", sep='\t', index_col=0, header=0)

    # Calculate variance for each gene (row)
    variances = df.var(axis=1)

    # Get the 5000 genes with the highest variance
    top_genes = df.loc[variances.nlargest(num_genes).index]

    # Save the filtered data
    top_genes.to_csv(f"results_assignment3/top_{num_genes}_variable_genes.tsv", sep='\t')
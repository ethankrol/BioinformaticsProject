import pandas as pd
from sklearn.cluster import AffinityPropagation
from sklearn.preprocessing import StandardScaler
import os

print("HI")

# Define directories
results_dir = "results_assignment3"
gene_data_dir = os.path.join(results_dir, "top_variable_genes")  # Assuming gene data files are stored here
metadata_path = os.path.join(results_dir, "filtered_metadata.tsv")

# Create results directory if it doesn't exist
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Define parameters
num_genes_list = [10, 100, 1000, 5000, 10000]
# Affinity Propagation does not require specifying the number of clusters, but you can adjust the 'preference' parameter
# to influence the number of clusters. Here, we keep it simple and let the algorithm determine it.

# Load metadata once
metadata = pd.read_csv(metadata_path, sep='\t', index_col=0, header=0)

# Iterate over different numbers of top variable genes
for num_genes in num_genes_list:
    # Construct the file path for the top variable genes
    gene_data_path = os.path.join(gene_data_dir, f"top_{num_genes}_variable_genes.tsv")
    
    # Check if the gene data file exists
    if not os.path.isfile(gene_data_path):
        print(f"Gene data file for top {num_genes} genes not found at {gene_data_path}. Skipping.")
        continue
    
    # Load gene expression data
    df = pd.read_csv(gene_data_path, sep='\t', index_col=0, header=0)
    
    # Transpose the data so that each row represents a sample
    df = df.T
    
    # Scale the data (mean=0, variance=1)
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(df)
    
    # Perform Affinity Propagation clustering
    # Note: AffinityPropagation may produce a large number of clusters by default.
    # You can adjust the 'preference' parameter to control the number of clusters.
    # A higher preference leads to more clusters, and vice versa.
    affinity = 'rbf'  # Other options: 'nearest_neighbors', 'precomputed'
    preference = -50  # Adjust this value based on your data. Lower values may result in fewer clusters.
    
    print(f"Performing Affinity Propagation clustering for top {num_genes} genes...")
    affinity_propagation = AffinityPropagation(affinity=affinity, preference=preference, random_state=1)
    clusters = affinity_propagation.fit_predict(df_scaled)
    
    # Number of clusters found
    num_clusters = len(set(clusters))
    print(f"Number of clusters found: {num_clusters}")
    
    # Prepare cluster assignments DataFrame
    clusters_df = pd.DataFrame({
        'refinebio_accession_code': df.index,
        'Cluster': clusters
    })
    
    # Merge cluster assignments with metadata
    metadata_with_clusters = metadata.merge(clusters_df, left_on='refinebio_accession_code', right_on='refinebio_accession_code')
    
    # Define output file path
    output_file = os.path.join(results_dir, f"affinity_propagation_clusters_from_{num_genes}_genes.tsv")
    
    # Save the clustered metadata to a TSV file
    metadata_with_clusters.to_csv(output_file, sep='\t', index=False)
    
    print(f"Clustered metadata saved to {output_file}\n")

print("Affinity Propagation clustering completed for all specified gene sets.")

import pandas as pd
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import StandardScaler

results_dir = "results_assignment3"

num_genes_list = [10, 100, 1000, 5000, 10000]
num_clusters_list = [2]

for num_genes in num_genes_list:
    for num_clusters in num_clusters_list:
        # Load gene data, transpose it for spectral clustering
        df = pd.read_csv(f"{results_dir}/top_{num_genes}_variable_genes.tsv", sep='\t', index_col=0, header=0)

        metadata = pd.read_csv(f"results/filtered_metadata.tsv", sep='\t', index_col=0, header=0)

        # transpose the data so that each row is a sample
        df = df.T

        # scale the data with mean 0 and variance 1
        scaler = StandardScaler()
        df_scaled = scaler.fit_transform(df) 

        # Perform spectral clustering
        spectral = SpectralClustering(n_clusters=num_clusters, affinity='nearest_neighbors', random_state=1)
        clusters = spectral.fit_predict(df_scaled)

        # add cluster info to metadata
        clusters_df = pd.DataFrame({'refinebio_accession_code': df.index, 'Cluster': clusters})
        metadata_with_clusters = metadata.merge(clusters_df, on='refinebio_accession_code')

        # write out results
        metadata_with_clusters.to_csv(f"{results_dir}/{num_clusters}_spectral_clusters_from_{num_genes}_genes.tsv", sep='\t', index=False)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score, classification_report

gene_10 = pd.read_csv(f"results_assignment3/K_means/top_10_variable_genes.tsv", sep='\t', index_col=0, header=0)
gene_100 = pd.read_csv(f"results_assignment3/K_means/top_100_variable_genes.tsv", sep='\t', index_col=0, header=0)
gene_1000 = pd.read_csv(f"results_assignment3/K_means/top_1000_variable_genes.tsv", sep='\t', index_col=0, header=0)
gene_5000 = pd.read_csv(f"results_assignment3/K_means/top_5000_variable_genes.tsv", sep='\t', index_col=0, header=0)
gene_10000 = pd.read_csv(f"results_assignment3/K_means/top_10000_variable_genes.tsv", sep='\t', index_col=0, header=0)
filtered_metadata = pd.read_csv(f"results/filtered_metadata.tsv", sep='\t', index_col=0, header=0)
kMean = pd.read_csv(f"results_assignment3/K_means/K_means_6_clusters_from_5000_genes.tsv", sep='\t', index_col=0, header=0)
timeStat = filtered_metadata['time_status']
kMean_Clusters = kMean['Cluster']

# Split into Training and Testing
X_train, X_test, y_train, y_test = train_test_split(gene_5000.T, timeStat, test_size=0.2, random_state=42)
# KNN Classifer
knn = KNeighborsClassifier(n_neighbors=10)
# Train the Model
knn.fit(X_train, y_train)
# Evaluate
y_pred = knn.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)
print(classification_report(y_test, y_pred))
pred = knn.predict(gene_5000.T)
cluster_predictions_df = pd.DataFrame(pred, index=timeStat.index, columns=['KNN_Label'])
cluster_predictions_df = cluster_predictions_df.rename_axis("sample")

# result = pd.read_csv(f"assignment_4/predictive_model_labels_5000_genes.tsv", sep='\t', index_col=0, header=0)
# result = pd.merge(result, cluster_predictions_df, on='sample', how='left')
# print(result)
# result.to_csv(f"assignment_4/predictive_model_labels_5000_genes.tsv", sep='\t')

# Split into Training and Testing
X_train, X_test, y_train, y_test = train_test_split(gene_5000.T, kMean_Clusters, test_size=0.2, random_state=42)
# KNN Classifer
knn = KNeighborsClassifier(n_neighbors=10)
# Train the Model
knn.fit(X_train, y_train)
# Evaluate
y_pred = knn.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)
print(classification_report(y_test, y_pred))
pred = knn.predict(gene_5000.T)
cluster_predictions_df = pd.DataFrame(pred, index=kMean_Clusters.index, columns=['Predicted Cluster'])


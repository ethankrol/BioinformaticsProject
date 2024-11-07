import pandas as pd
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.preprocessing import StandardScaler

# Load the expression data
data = pd.read_csv('results_assignment3/top_5000_variable_genes.tsv', sep='\t')
data = data.transpose()

# Load and merge the label data with expression data
label_data = pd.read_csv("results_assignment3/GMM/GMM_5000_genes.csv")[["refinebio_accession_code", "time_status", "cluster"]]
data = pd.merge(data, label_data, left_index=True, right_on='refinebio_accession_code')

# Only process the 'time_status' group
group = "time_status"

# Convert "day" to 0 and "month" to 1 in the "time_status" group
data[group] = data[group].replace({"day": 0, "month": 1})

# Define gene subset sizes
gene_subset_sizes = [10, 100, 1000, 10000]

for n_genes in gene_subset_sizes:
    # Select the first `n_genes` features
    X = data.iloc[:, 1:n_genes + 1]  # Adjust for index offset
    X.columns = X.columns.astype(str)  # Ensure columns are strings

    # Define the target labels (y) for the 'time_status' group
    y = data[group]

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

    # Standardize features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Initialize and train the Naive Bayes classifier
    nb = GaussianNB()
    nb.fit(X_train, y_train)

    # Make predictions on the test set
    y_pred = nb.predict(X_test)

    # Print classification metrics
    print(f"\nResults for Naive Bayes classifier on group: {group} with {n_genes} genes")
    print(classification_report(y_test, y_pred))

    # Calculate AUC for binary classification
    y_proba = nb.predict_proba(X_test)[:, 1]  # Probability for positive class
    auc_score = roc_auc_score(y_test, y_proba)
    print(f"AUC Score for {group} with {n_genes} genes: {auc_score:.4f}")

    # Write out model labels for all samples in 'time_status'
    model_labels = nb.predict(X)
    model_output = pd.DataFrame({"sample": data["refinebio_accession_code"], f"{group}_NB_label_{n_genes}_genes": model_labels})
    model_output.to_csv(f"assignment_4/naive_bayes_{n_genes}_genes_time_status.tsv", sep='\t', index=False)

    # Step 2e: Predict clusters with Naive Bayes
    y_cluster = data["cluster"]

    # Split the data for cluster prediction
    X_train, X_test, y_train, y_test = train_test_split(X, y_cluster, test_size=0.2, random_state=1)

    # Standardize and train Naive Bayes on the clusters
    nb.fit(X_train, y_train)
    y_pred = nb.predict(X_test)

    # Print classification metrics for cluster prediction
    print(f"\nResults for Naive Bayes classifier on group: 'cluster' with {n_genes} genes")
    print(classification_report(y_test, y_pred))

    # Write out model labels for all samples in 'cluster'
    model_labels_cluster = nb.predict(X)
    model_output_cluster = pd.DataFrame({"sample": data["refinebio_accession_code"], f"cluster_NB_label_{n_genes}_genes": model_labels_cluster})
    model_output_cluster.to_csv(f"assignment_4/Naive_Bayes/naive_bayes_{n_genes}_genes_cluster.tsv", sep='\t', index=False)

print("Naive Bayes model results saved for 'time_status' and 'cluster' predictions at all gene levels.")

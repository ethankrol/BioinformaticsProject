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

# Prepare the features (X), excluding the last three columns with label information
X = data.iloc[:, 1:-3]
X.columns = X.columns.astype(str)  # Ensure columns are strings

# Only process the 'time_status' group
group = "time_status"

# Define the target labels (y) for the 'time_status' group
y_time_status = data[group]

# Convert "day" to 0 and "month" to 1 in the "time_status" group
y_time_status = y_time_status.replace({"day": 0, "month": 1})

# Split the data into training and testing sets for 'time_status'
X_train, X_test, y_train, y_test = train_test_split(X, y_time_status, test_size=0.2, random_state=1)

# Standardize features
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Initialize and train the Naive Bayes classifier for 'time_status'
nb = GaussianNB()
nb.fit(X_train, y_train)

# Make predictions on the test set
y_pred = nb.predict(X_test)

# Print classification metrics
print(f"Results for Naive Bayes classifier on group: {group}")
print(classification_report(y_test, y_pred))

# Calculate AUC for binary classification
y_proba = nb.predict_proba(X_test)[:, 1]  # Probability for positive class
auc_score = roc_auc_score(y_test, y_proba)
print(f"AUC Score for {group}: {auc_score}")

# Write out model labels for all samples in 'time_status'
model_labels_time_status = nb.predict(X)
model_output_time_status = pd.DataFrame({"sample": data["refinebio_accession_code"], f"{group}_NB_label": model_labels_time_status})
model_output_time_status.to_csv(f"assignment_4/naive_bayes_5000_genes_time_status.tsv", sep='\t', index=False)

# Step 2e: Predict clusters with Naive Bayes

# Define the target labels (y) for the 'cluster' group
y_cluster = data["cluster"]

# Split the data for cluster prediction
X_train, X_test, y_train, y_test = train_test_split(X, y_cluster, test_size=0.2, random_state=1)

# Standardize and train Naive Bayes on the clusters
nb.fit(X_train, y_train)
y_pred = nb.predict(X_test)

# Print classification metrics for cluster prediction
print(f"Results for Naive Bayes classifier on group: 'cluster'")
print(classification_report(y_test, y_pred))

# Write out model labels for all samples in 'cluster'
model_labels_cluster = nb.predict(X)
model_output_cluster = pd.DataFrame({"sample": data["refinebio_accession_code"], "cluster_NB_label": model_labels_cluster})
model_output_cluster.to_csv(f"assignment_4/naive_bayes_5000_genes_cluster.tsv", sep='\t', index=False)

print("Naive Bayes model results saved for both 'time_status' and 'cluster' predictions.")

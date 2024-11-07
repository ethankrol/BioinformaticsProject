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
y = data[group]

# Convert "day" to 0 and "month" to 1 in the "time_status" group
y = y.replace({"day": 0, "month": 1})

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
print(f"Results for Naive Bayes classifier on group: {group}")
print(classification_report(y_test, y_pred))

# Calculate AUC for binary classification
y_proba = nb.predict_proba(X_test)[:, 1]  # Probability for positive class
auc_score = roc_auc_score(y_test, y_proba)
print(f"AUC Score for {group}: {auc_score}")

# Write out model labels for all samples
model_labels = nb.predict(X)
model_output = pd.DataFrame({"sample": data["refinebio_accession_code"], f"{group}_NB_label": model_labels})
model_output.to_csv(f"assignment_4/naive_bayes_5000_genes.tsv", sep='\t', index=False)

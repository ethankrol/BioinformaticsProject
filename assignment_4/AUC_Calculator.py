import pandas as pd
from sklearn.metrics import roc_auc_score

# Load your data
df = pd.read_csv('assignment_4/predictive_model_labels_5000_genes.tsv', sep='\t') # Uncomment if loading from a CSV

# Example data setup: Replace this with your data


# Convert categorical labels to binary (0 and 1)
label_mapping = {'One_Day_After_Defeat': 0, 'One_Month_After_Defeat': 1}
for col in df.columns:
    df[col] = df[col].map(label_mapping)

# Separate true labels and model predictions
y_true = df['time_status_NB_label']
model_columns = ['SVM_label', 'RFC_label', 'LogisticRegression_label', 'KNN_label', 'NaiveBayes_label']

# Calculate AUC for each model
auc_scores = {}
for model in model_columns:
    auc = roc_auc_score(y_true, df[model])
    auc_scores[model] = auc
    print(f"AUC for {model}: {auc}")

print("AUC Scores:", auc_scores)

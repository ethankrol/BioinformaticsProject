import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score, roc_curve, auc
from sklearn.preprocessing import StandardScaler

data = pd.read_csv('results_assignment3/top_5000_variable_genes.tsv', sep = '\t')  

data = data.transpose()

label_data = pd.read_csv("results_assignment3/GMM/GMM_5000_genes.csv")[["refinebio_accession_code", "time_status", "cluster"]]

data = pd.merge(data, label_data, left_index=True, right_on='refinebio_accession_code')

X = data.iloc[:, 1:-3]
scaler = StandardScaler()
X = scaler.fit_transform(X)

num_genes = 5000

for group in ["time_status"]:

    y = data[group]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=3)

    rfc = RandomForestClassifier()
    rfc.fit(X_train, y_train)

    y_pred = rfc.predict(X_test)
    print(f"Group: {group}")
    print(classification_report(y_test, y_pred))   


model_labels = rfc.predict(X)
print(model_labels)
model_output = pd.read_csv("assignment_4/predictive_model_labels_5000_genes.tsv", sep = '\t')
model_output["RFC_label"] = model_labels
model_output.to_csv(f"assignment_4/predictive_model_labels_{num_genes}_genes.tsv", sep="\t", index=False)
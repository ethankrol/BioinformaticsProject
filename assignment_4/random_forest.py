import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score, roc_curve, auc
from sklearn.preprocessing import StandardScaler

num_genes_list = [10, 100, 1000, 5000, 10000]

all_results = []

for num_genes in num_genes_list:
    data = pd.read_csv(f'results_assignment3/top_{num_genes}_variable_genes.tsv', sep = '\t')  

    data = data.transpose()

    label_data = pd.read_csv(f"results_assignment3/GMM/GMM_{num_genes}_genes.csv")[["refinebio_accession_code", "time_status", "cluster"]]

    data = pd.merge(data, label_data, left_index=True, right_on='refinebio_accession_code')

    X = data.iloc[:, 1:-3]

    group = "time_status"
    y = data[group]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    rfc = RandomForestClassifier(random_state = 1)
    rfc.fit(X_train, y_train)

    y_pred = rfc.predict(X_test)
    y_prob = rfc.predict_proba(X_test)[:, 1]
    print(f"Group: {num_genes}")
    report = classification_report(y_test, y_pred, output_dict=True)
    print(report)  

    auc_score = roc_auc_score(y_test, y_prob) 

    metrics = {
        'num_genes': num_genes,
        'AUC': auc_score,
        'precision': report['weighted avg']['precision'],
        'recall': report['weighted avg']['recall'],
        'f1-score': report['weighted avg']['f1-score'],
        'support': report['weighted avg']['support'],
    }

    all_results.append(metrics)

    model_labels = rfc.predict(X)

    data[f'RFC_label_{num_genes}_genes'] = model_labels
    data[['refinebio_accession_code', f'RFC_label_{num_genes}_genes']].to_csv(
        f"assignment_4/predictive_model_labels_{num_genes}_genes.tsv", sep='\t', index=False
    )

results_df = pd.DataFrame(all_results)

results_df.to_csv('assignment_4/random_forest_results_all_groups.tsv', sep='\t', index=False)

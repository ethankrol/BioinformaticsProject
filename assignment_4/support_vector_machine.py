import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import StandardScaler

num_genes_list = [10, 100, 1000, 10000]

for num_genes in num_genes_list:
    # Load the data
    data = pd.read_csv(f"results_assignment3/top_{num_genes}_variable_genes.tsv", sep="\t")

    # Transpose to have samples as rows and genes as columns
    data = data.transpose()

    # Sample labels
    label_data = pd.read_csv(f"results_assignment3/2_spectral_clusters_from_{num_genes}_genes.tsv", sep="\t")[["refinebio_accession_code", "time_status", "Cluster"]]

    # Merge the data
    data = pd.merge(data, label_data, left_index=True, right_on="refinebio_accession_code")

    # create X data
    X = data.iloc[:, :-3]
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    for group in ["time_status"]:
        # create y data
        y = data[group]

        # split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

        # train svm
        svm = SVC()
        svm.fit(X_train, y_train)

        # predict/evaluate
        y_pred = svm.predict(X_test)
        print(classification_report(y_test, y_pred))

# write out model labels for each sample
model_labels = svm.predict(X)
print(model_labels)
model_output = pd.DataFrame(data={"sample": data["refinebio_accession_code"], "SVM_label": model_labels})
model_output.to_csv(f"assignment_4/predictive_model_labels_{num_genes}_genes.tsv", sep="\t", index=False)
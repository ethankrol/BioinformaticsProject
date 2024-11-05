import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import StandardScaler

num_genes = 5000

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

for group in ["time_status", "Cluster"]:
    # create y data
    y = data[group]

    # split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

    # standardize features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # train svm
    svm = SVC()
    svm.fit(X_train, y_train)

    # predict/evaluate
    y_pred = svm.predict(X_test)
    print(f"Group: {group}")
    #print("accuracy: ", accuracy_score(y_test, y_pred))
    print(classification_report(y_test, y_pred))
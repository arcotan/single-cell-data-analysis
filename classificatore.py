import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
import scanpy as sc
import numpy as np

MARKERS_PATH = "./results/markers_10X_P7_4_seurat.csv"
LABELS_PATH = "./filtered_dataset/tabulamuris/labels_10X_P7_4.csv"
markers_df = pd.read_csv(MARKERS_PATH)  
top_genes = markers_df.gene.unique()

adata = sc.read_10x_mtx(
    './filtered_dataset/tabulamuris/data_10X_P7_4/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=False
)

X = adata[:, top_genes].X.toarray()
labels_df = pd.read_csv(LABELS_PATH)
# y = labels_df['cluster.ids'].to_numpy()

# join between adata.obj_names and labels_df.cell
y = np.array([labels_df[labels_df.cell == cell]['cluster.ids'] for cell in adata.obs_names])


a = np.isnan(y).reshape(-1)
# remove the nan from both X that is a numpy array
y = y[a == False]
X = X[~a, :]

 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, stratify=y)
clf = RandomForestClassifier()
clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)
print(f1_score(y_test, y_pred, average='macro'))

importance = clf.feature_importances_

a = 1
# Read the data
#data = pd.read_csv('./scanpy/.csv')
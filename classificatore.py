import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score
import scanpy as sc
import numpy as np

#MARKERS_PATH = "./results/markers_10X_P7_4_seurat.csv"
MARKERS_PATH = './results/draft_rank_genes_scanpy.csv'

LABELS_PATH = "./filtered_dataset/tabulamuris/labels_10X_P7_4.csv"
markers_df = pd.read_csv(MARKERS_PATH)  
top_genes = markers_df.gene.unique()

adata = sc.read_10x_mtx(
    './filtered_dataset/tabulamuris/data_10X_P7_4/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=False
)

#X = adata[:, top_genes].X.toarray()
X=adata.X.toarray()
labels_df = pd.read_csv(LABELS_PATH)
# y = labels_df['cluster.ids'].to_numpy()

# join between adata.obj_names and labels_df.cell
y = np.array([labels_df[labels_df.cell == cell]['cluster.ids'] for cell in adata.obs_names])


a = np.isnan(y).reshape(-1)
# remove the nan from both X that is a numpy array
y = y[a == False]
X = X[~a, :]

from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_predict

clf = RandomForestClassifier()
y_pred = cross_val_predict(clf, X, y, cv=5)
report = classification_report(y, y_pred, output_dict=True)
clf.fit(X, y)
feature_importance = clf.feature_importances_ # calcolare importances di ciascun fold e aggregare i rankings?

all_features_list = adata.var_names

# sort both lists by feature_importance mantaining the same order
sorted_list = sorted(zip(feature_importance, all_features_list), reverse=True)

# top 120 genes with feature importance different from all the others

current_score = sorted_list[120][0]
top_sorted_list = [x[1] for x in sorted_list if x[0] >= current_score]
intersection = set(top_genes).intersection(set(top_sorted_list))



# all genes
# {'precision': 0.942666699684296, 'recall': 0.9415041782729805, 'f1-score': 0.9404974321493618, 'support': 359}
# top
# {'precision': 0.9700833699346517, 'recall': 0.9693593314763231, 'f1-score': 0.9693161075050654, 'support': 359}
print('fine')
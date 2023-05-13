import argparse
import pandas as pd
import scanpy as sc
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_predict

"""
parser = argparse.ArgumentParser(
    description='Markers validation through classification'
)
parser.add_argument(
    '--markers-path',
    action='store',
    type=str,
    required=True,
    help="The path to the .csv file with markers."
)
parser.add_argument(
    '--data-path',
    action='store',
    type=str,
    required=True,
    help="The path to the dataset."
)
parser.add_argument(
    '--labels-path',
    action='store',
    type=str,
    required=True,
    help="The path to the .csv file with labels."
)
parser.add_argument(
    '--out-path',
    action='store',
    type=str,
    required=True,
    help="The path to the output folder."
)
parser.add_argument(
    '--simulated',
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Wheter the dataset is simulated."
)
args = parser.parse_args()
data_path = args.data_path
markers_path = args.markers_path
labels_path = args.labels_path
out_path = args.out_path
simulated = args.simulated
"""

# -------- To test it without cmd line args --------
# TODO: change paths
markers_path = "./results/seurat/markers_10X_P7_4_seurat.csv"
data_path = "./filtered_dataset/tabulamuris/data_10X_P7_4/"
labels_path = "./filtered_dataset/tabulamuris/labels_10X_P7_4.csv"
out_path = "./"
simulated = False
# --------------------------------------------------

def apply_classifier(X, y):
    clf = RandomForestClassifier()
    y_pred = cross_val_predict(clf, X, y, cv=5)
    report = classification_report(y, y_pred, output_dict=True)
    clf.fit(X, y)
    return report, clf.feature_importances_

if simulated:
    pass
    # TODO:
    # read matrix
    # adata = ad.AnnData(matrix)
    # adata.obs_names = ...
    # adata.var_names = ...
else:
    adata = sc.read_10x_mtx(
        data_path,
        var_names='gene_symbols',
        cache=False
    )

markers_df = pd.read_csv(markers_path)
markers = markers_df.gene.unique()

y_df = pd.read_csv(labels_path)
y = np.array([y_df[y_df.cell == cell]['cluster.ids'] for cell in adata.obs_names])
nan_mask = ~np.isnan(y).reshape(-1)
y = y[nan_mask]

X_all = adata.X.toarray()[nan_mask, :]
X_top = adata[:, markers].X.toarray()[nan_mask, :]

report_all, feature_importance = apply_classifier(X_all, y)
report_top, _ = apply_classifier(X_top, y)

pd.DataFrame(report_all).transpose().to_csv(out_path+"clf_score_all.csv")
pd.DataFrame(report_top).transpose().to_csv(out_path+"clf_score_top.csv")

sorted_idx = (-feature_importance).argsort()
features_sorted = adata.var_names[sorted_idx]
importaces_sorted = feature_importance[sorted_idx]
pd.DataFrame(
    {'genes' : features_sorted, 'importaces' : importaces_sorted}
    ).to_csv(out_path+"importances.csv")

# TODO:
# - valutare intersezione
# - valutare bontà del ranking allenando con markers più in basso nella classifica?
# - selezionare il numero di marcatori in base alle prestazioni del classificatore
"""
# select top features with recursive feature elimination and random forest
from sklearn.feature_selection import RFE
selector = RFE(RandomForestClassifier(), n_features_to_select=120, step=0.5)
selector.fit(X_all, y)
important_features = adata.var_names[selector.support_] # to test!
"""

important_features = features_sorted[0:120]
# important_features = [f for f, i in zip(features_sorted, importaces_sorted) if i >= importaces_sorted[119]]
intersection = set(markers).intersection(set(important_features))

#        rank di randomforest                  rank del tool
# gene1         100 *                         * 1 - (0-20)
# gene2         1 *                            
# gene3         1 *                            
# gene4         0 *                            
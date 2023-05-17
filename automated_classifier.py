import pandas as pd
import numpy as np
import scanpy as sc
import os
import sys
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.metrics import f1_score
from sklearn.model_selection import cross_val_predict
import pickle


DATASET_TAGS = ['tabula-muris-heart']#, 'tabula-muris-marrow_P7_3', 'peripheal-blood', 'kumar-4-hard', 'kumar-8-hard']
N_MARKERS = 50#50

def apply_classifier(X, y):
	clf = RandomForestClassifier()
	y_pred = cross_val_predict(clf, X, y, cv=5)
	f1 = f1_score(y, y_pred)
	return f1

# assuming cells in y_df >= cells in adata and cluster.ids != nan
def read_data(data_path, labels_path, markers_path):
	adata = sc.read_10x_mtx(
		data_path,
		var_names='gene_symbols',
		cache=False
	)
	y_df = pd.read_csv(labels_path, index_col=0)
	y_df = pd.DataFrame(adata.obs_names, columns=["cell"]).join(y_df, on="cell")
	y = np.array(y_df['cluster.ids'])
	markers_df = pd.read_csv(markers_path)
	return adata, y, markers_df

def process(tool, clusters, weights, markers_df, adata, y_ovr):
    scores = {}
    m = [] # [[f1_0, ...., f1_50]_1, ...,[f1_0, ...., f1_50]_n_clusters]
    for cluster in clusters:
        for i in range(1, N_MARKERS+1):
            markers = markers_df[
                (markers_df['cluster']==cluster) & 
                (markers_df['tool']==tool) & 
                (markers_df['rank']<=i)
                ]['gene'].values
            X = adata[:, markers].X.toarray()
            f1 = apply_classifier(X, y_ovr[cluster])
            if tool not in scores:
                scores[tool] = {}
            if cluster not in scores[tool]:
                scores[tool][cluster] = []
            scores[tool][cluster].append(f1)
        
        m.append(scores[tool][cluster])
        # mean for each cluster
        print("Cluster {} done for tool {}".format(cluster, tool), file=sys.stderr)
    
    m = np.array(m)

    # multiply each column of m pointwise with the weights array
    
    assert m.shape[0] == len(weights)
    assert m.shape[1] == N_MARKERS

    m = np.multiply(m, weights.reshape(-1, 1))

    scores[tool]['mean'] = np.mean(m, axis=0) / np.sum(weights)
        
    return scores

# scores = {
# 	'monocle': {    # una sola chiave
# 		'1': [f1_0, ..., f1_50],
#       '2' : [f1_0, ..., f1_50],
#       'mean': [meanf1_0, ..., meanf1_50]
# 	}
# }


# ('monocle', {})

for dataset in DATASET_TAGS:
	data_path = "./dataset/{}-filtered/10X/".format(dataset)
	labels_path = "./dataset/{}-filtered/labels.csv".format(dataset)
	markers_path = "./results/aggregate/{}/markers.csv".format(dataset)
	out_path = "./results/aggregate/{}/".format(dataset)
	scores_path = './results/{}/clf_scores.pickle'.format(dataset)

	adata, y, markers_df = read_data(data_path, labels_path, markers_path)
	clusters, weights = np.unique(y, return_counts=True)
	n_clusters = len(clusters)
	tools = markers_df['tool'].unique()

	y_ovr = {}
	for cluster in clusters:
		y_ovr[cluster] = np.array(y==cluster, dtype=int)

	# train binary classifiers for each cluster using tools markers and compute f1
	f1_markers = {}
	for tool in tools:
		f1_markers_tool = []
		for cluster in clusters:
			markers = markers_df[(markers_df['cluster']==cluster) & (markers_df['tool']==tool)]['gene'].values
			X_markers = adata[:, markers].X.toarray()
			f1 = apply_classifier(X_markers, y_ovr[cluster])
			f1_markers_tool.append(f1)
		f1_markers[tool] = round((weights*np.array(f1_markers_tool)).sum()/weights.sum(), 5)
	
	# train binary classifiers for each cluster using all genes and compute f1
	f1_all = []
	X_all = adata.X.toarray()
	for cluster in clusters:
		f1 = apply_classifier(X_all, y_ovr[cluster])
		f1_all.append(f1)
	f1_markers['all'] = round((weights*np.array(f1_all)).sum()/weights.sum(), 5)
	pd.DataFrame.from_dict(
		data=f1_markers,
		orient='index').to_csv(out_path+"clf_f1.csv", header=['f1'])
	
	# rank genes using recursive feature elimination
	selector = RFE(RandomForestClassifier(), n_features_to_select=N_MARKERS, step=0.5)
	cluster_features = {}
	for cluster in clusters:
		selector.fit(X_all, y_ovr[cluster])
		sorted_idx = (selector.ranking_).argsort()
		rfe_features_sorted = adata.var_names[sorted_idx]
		cluster_features[cluster] = rfe_features_sorted
	pd.DataFrame(cluster_features).to_csv(out_path+"rfe_ranking.csv")

	if not os.path.exists(scores_path):
		scores = Parallel(n_jobs=len(tools))(delayed(process)( # TODO: vedere se funziona
			tool,
			clusters,
			weights,
			markers_df,
			adata,
			y_ovr) for tool in tools)
		scores2 = {}
		for i, el in enumerate(scores):
			scores2[list(el.keys())[0]] = el[list(el.keys())[0]]
		scores = scores2
		with open(scores_path, 'wb') as handle:
			pickle.dump(scores, handle, protocol=pickle.HIGHEST_PROTOCOL)   

# plot with increasing features (feature di un classificatore in particolare errata...)
# multi-class?

# for cluster in clusters:
# 	rfe = RFECV(estimator=RandomForestClassifier())
# 	rfe.fit(X_all, y_ovr[cluster])
# 	rfe.show()
import pandas as pd
import numpy as np
import scanpy as sc
import sys
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.metrics import f1_score
from sklearn.model_selection import cross_val_predict
import pickle

# DATASET_TAGS = ['tabula-muris-heart', 'tabula-muris-marrow_P7_3', 'peripheral-blood', 'zheng-4', 'zheng-8']
DATASET_TAGS = ['PBMC1', 'PBMC2']
N_MARKERS = 500
NULL_MODEL_ATTEMPTS = 50

def apply_classifier(X, y):
	clf = RandomForestClassifier(n_jobs=4)
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
	y = np.array(y_df['cluster.ids']).astype(int)
	markers_df = pd.read_csv(markers_path)
	return adata, y, markers_df

def process(tool, clusters, weights, markers_df, adata, y_ovr):
    scores = {}
    f1_lists = []
    for cluster in clusters:
        for i in range(1, N_MARKERS+1):
            markers = markers_df[
                (markers_df['cluster']==cluster) & 
                (markers_df['tool']==tool) & 
                (markers_df['rank']<=i)
                ]['gene'].values
            X = adata[:, markers].X.toarray()
            f1 = apply_classifier(X, y_ovr[cluster])
            if cluster not in scores:
                scores[cluster] = []
            scores[cluster].append(f1)
        f1_lists.append(scores[cluster])
        print("Cluster {} done for tool {}".format(cluster, tool), file=sys.stderr)
    scores['mean'] = np.sum(np.array(f1_lists)*weights.reshape(weights.size, 1), axis=0)/np.sum(weights)
    return (tool, scores)

#TODO
def weighted_std(values, weights):
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return np.sqrt(variance)

for dataset in DATASET_TAGS:
	print("# Processing dataset {}".format(dataset))
	data_path = "./dataset/{}-Filtered/10X/".format(dataset)
	labels_path = "./dataset/{}-Filtered/celltypist_labels.csv".format(dataset)
	markers_path = "./results/aggregate_new_pbmcs/{}/markers.csv".format(dataset)
	out_path = "./results/aggregate_new_pbmcs/{}/".format(dataset)
	scores_path = './results/{}/clf_scores.pickle'.format(dataset)

	adata, y, markers_df = read_data(data_path, labels_path, markers_path)
	na_id = np.isnan(y)
	adata = adata[~na_id]
	y = y[~na_id]

	clusters, weights = np.unique(y, return_counts=True)
	n_clusters = len(clusters)
	tools = markers_df['tool'].unique()

	y_ovr = {}
	for cluster in clusters:
		y_ovr[cluster] = np.array(y==cluster, dtype=int)

	# train binary classifiers for each cluster using tools markers and compute f1
	print("## Training 1 vs Rest RF on markers")
	f1_std = {}
	f1_markers = {}
	for tool in tools:
		print("Tool {}".format(tool))
		f1_markers_tool = []
		for cluster in clusters:
			markers = markers_df[(markers_df['cluster']==cluster) & (markers_df['tool']==tool)]['gene'].values
			X_markers = adata[:, markers].X.toarray()
			f1 = apply_classifier(X_markers, y_ovr[cluster])
			f1_markers_tool.append(f1)
			print("Cluster {} done".format(cluster))
		f1_markers[tool] = round((weights*np.array(f1_markers_tool)).sum()/weights.sum(), 5)
		f1_std[tool] = weighted_std(f1_markers_tool, weights)
	
	# train binary classifiers for each cluster using all genes and compute f1
	print("## Training 1 vs Rest RF on all features")
	f1_all = []
	X_all = adata.X.toarray()
	for cluster in clusters:
		f1 = apply_classifier(X_all, y_ovr[cluster])
		f1_all.append(f1)
		print("Cluster {} done".format(cluster))
	f1_markers['all'] = round((weights*np.array(f1_all)).sum()/weights.sum(), 5)
	f1_std['all'] = weighted_std(f1_all, weights)
	pd.DataFrame.from_dict(
		data=f1_markers,
		orient='index').to_csv(out_path+"clf_f1.csv", header=['f1'])
	
	# train binary classifiers for each cluster using all genes and compute f1
	print(f"## Training 1 vs Rest RF on a random subset of {N_MARKERS} features")
	f1_sample_std = []
	f1_sample = []
	for i in range(NULL_MODEL_ATTEMPTS):
		f1_sample_attempt = []
		X_sample = adata[:, np.random.choice(adata.shape[1], N_MARKERS, replace=False)].X.toarray()
		for cluster in clusters:
			f1 = apply_classifier(X_sample, y_ovr[cluster])
			f1_sample_attempt.append(f1)
			print("Cluster {} done".format(cluster))
		f1_sample_std.append(weighted_std(f1_sample_attempt, weights))
		f1_sample.append(round((weights*np.array(f1_sample_attempt)).sum()/weights.sum(), 5))
	f1_markers['sample'] = round(np.mean(f1_sample), 5)
	f1_std['sample'] = round(np.mean(f1_sample_std), 5)
	pd.DataFrame.from_dict(
		data=f1_markers,
		orient='index').to_csv(out_path+"clf_f1_sample.csv", header=['f1'])
	pd.DataFrame.from_dict(
		data=f1_std,
		orient='index').to_csv(out_path+"clf_f1_std.csv", header=['f1'])
	
	# rank genes using recursive feature elimination
	print("## Ranking genes with RFE and RF")
	selector = RFE(RandomForestClassifier(n_jobs=8), n_features_to_select=N_MARKERS, step=0.5)
	cluster_features = {}
	for cluster in clusters:
		selector.fit(X_all, y_ovr[cluster])
		sorted_idx = (selector.ranking_).argsort()
		rfe_features_sorted = adata.var_names[sorted_idx]
		cluster_features[cluster] = rfe_features_sorted
		print("Cluster {} done".format(cluster))
	pd.DataFrame(cluster_features).to_csv(out_path+"rfe_ranking.csv")

	print("## Training 1 vs Rest RF with increasing # of features")
	scores = Parallel(n_jobs=len(tools))(delayed(process)(tool, clusters, weights, markers_df, adata, y_ovr) for tool in tools)
	scores_to_write = {}
	for (tool, score) in scores:
		scores_to_write[tool] = score
	with open(scores_path, 'wb') as fp:
		pickle.dump(scores_to_write, fp, protocol=pickle.HIGHEST_PROTOCOL)
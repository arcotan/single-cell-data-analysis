import pandas as pd

markers_lists = {
	'endothelial cell': ['Fabp4', 'Cdh5', 'Cav1'],
	'fibroblast': ['Ddr2', 'Tcf21', 'Col3a1', 'Col1a2', 'Col1a1'],
	'cardiac muscle cell': ['Nppa', 'Myl7', 'Sln'],
	'endocardial cell': ['Npr3', 'Pecam1'],
	'immune_cells': ['C1qa', 'H2-Eb1'],
	'myofibroblasts / smooth muscle': ['Myh11', 'Tcf21']
}

df = pd.read_csv("./results/aggregate/tabula-muris-heart/markers.csv")

ranks = pd.DataFrame(columns=['gene', 'cluster', 'rank', 'tool', 'truth_label'])

for cluster, marker_list in markers_lists.items():
	for marker in marker_list:
		to_concat = df[df['gene']==marker].copy()
		to_concat['truth_label'] = cluster
		ranks = pd.concat([ranks, to_concat])

ranks.to_csv("./results/aggregate/tabula-muris-heart/truth_markers.csv", index=False)
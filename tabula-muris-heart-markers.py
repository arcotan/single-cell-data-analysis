import pandas as pd

markers_lists = {
	'endothelial cell': ['Fabp4', 'Cdh5', 'Cav1'],
	'fibroblast': ['Ddr2', 'Tcf21', 'Col3a1', 'Col1a2', 'Col1a1', 'Myh11', 'Tcf21'],
	'cardiac muscle cell': ['Nppa', 'Myl7', 'Sln'],
	'endocardial cell': ['Npr3', 'Pecam1'],
	'immune cell': ['C1qa', 'H2-Eb1'],
}

df = pd.read_csv("./results/aggregate/tabula-muris-heart/markers.csv")

ranks = pd.DataFrame(columns=['gene', 'cluster', 'rank', 'tool', 'truth_label'])

for cluster, marker_list in markers_lists.items():
	for marker in marker_list:
		to_concat = df[df['gene']==marker].copy()
		to_concat['truth_label'] = cluster
		ranks = pd.concat([ranks, to_concat])

ranks.to_csv("./results/aggregate/tabula-muris-heart/truth_markers.csv", index=False)


import pickle
import matplotlib.pyplot as plt
import numpy as np

TOOL_COLOR = {
    "COTAN": "#66c2a5", 
    "seurat": "#fc8d62", 
    "scanpy": "#8da0cb", 
    "monocle": "#e78ac3", 
    "scvi": "#a6d854"
}

TOOL_TAGS = {
	'COTAN': 'COTAN',
	'seurat': 'Seurat',
	'scanpy': 'Scanpy',
	'monocle':'Monocle',
	'scvi-tools': 'scVI-tools',
}

order = ['COTAN', 'Seurat', 'Scanpy', 'Monocle', 'scVI-tools']

DATASET_TAGS = ['tabula-muris-heart', 'tabula-muris-marrow_P7_3', 'peripheal-blood', 'zheng-4', 'zheng-8']
N_MARKERS = 50

for dataset in DATASET_TAGS:
	scores_path = './results/{}/clf_scores.pickle'.format(dataset)
	with open(scores_path, 'rb') as handle:
		scores = pickle.load(handle)
	tools = scores.keys()
	fig, ax = plt.subplots(2, 3, figsize=(21, 9))
	fig.suptitle(dataset.replace('-', ' ').title(), fontsize=16)
	min_val = np.min([np.min(scores[tool]['mean']) for tool in tools])
	max_val = np.max([np.max(scores[tool]['mean']) for tool in tools])
	for j, tool in enumerate(tools): 
		print(tool)
		i = 1-j%2 if j< 5 else 2			
		max_i = np.argmax(scores[tool]['mean'])
		ax[i%2, j//2].plot([i for i in range(1, N_MARKERS+1)], scores[tool]['mean'], marker='.')
		
		# print max max in red
		ax[i%2, j//2].plot([max_i+1], [scores[tool]['mean'][max_i]], marker='o', color='red')
		
		# annotation for max
		ax[i%2, j//2].annotate('max: {}'.format(round(scores[tool]['mean'][max_i], 3)), xy=(max_i+1, scores[tool]['mean'][max_i]), xytext=(max_i+1, scores[tool]['mean'][max_i]+0.01), arrowprops=dict(facecolor='black', shrink=0.05))
		
		# annotation for second max
		
		# Only last two plots have x label (it's the same for all)
		if (i>1):
			ax[i%2, j//2].set_ylabel("f1 weighted")
			ax[i%2, j//2].set_xlabel("# features")
		# print x tick for max 
		ax[i%2, j//2].plot([max_i+1, max_i+1], [min_val, scores[tool]['mean'][max_i]], linestyle='--', color='red')
		ax[i%2, j//2].set_xticks([max_i+1], [max_i+1], minor=True)
		# Same scale for all plots
		ax[i%2, j//2].set_yticks([round(i, 2) for i in np.arange(min_val, 1, (0.005 if min_val > 0.85 else 0.05))])
		# Set title
		ax[i%2, j//2].set_title(tool)
		ax[i%2, j//2].grid()
	ax[0, 2].set_visible(False)
	plt.savefig('./results/{}/clf_plots.eps'.format(dataset))
	plt.savefig('./results/{}/clf_plots.png'.format(dataset))

	
best_for_COTAN = (0, '')
worst_for_COTAN = (0, '')
for dataset in DATASET_TAGS:
    scores_path = './results/{}/clf_scores.pickle'.format(dataset)
    with open(scores_path, 'rb') as handle:
        scores = pickle.load(handle)
    tools = scores.keys()
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))
    fig.suptitle(dataset.replace('-', ' ').title(), fontsize=16)
    for j, tool in enumerate(tools):
        ax.plot([i for i in range(1, N_MARKERS+1)], scores[tool]['mean'], color=TOOL_COLOR[tool], label=TOOL_TAGS[tool])
        
        # dataset where COTAN on average is the best compared to other tools
        if best_for_COTAN == (0, '') or np.average(scores['COTAN']['mean']) > np.average(scores[tool]['mean']):
            best_for_COTAN = (np.average(scores['COTAN']['mean']), dataset)
        # dataset where COTAN on average is the worst compared to other tools
        if worst_for_COTAN == (0, '') or np.average(scores['COTAN']['mean']) < np.average(scores[tool]['mean']):
            worst_for_COTAN = (np.average(scores['COTAN']['mean']), dataset)
    ax.set_ylabel("f1 weighted")
    ax.set_xlabel("# features")
    # ax.set_yticks([round(i, 2) for i in np.arange(min_all[dataset], max_all[dataset], (0.005 if min_val > 0.9 else 0.05))])
    # ax.legend(fontsize='20')
    handles, labels = plt.gca().get_legend_handles_labels()
    ax.legend([handles[order.index(index)] for index in order],[labels[labels.index(index)] for index in order], fontsize='20')
    ax.grid()
    
    plt.savefig('./results/{}/clf_ALL.eps'.format(dataset))
    plt.savefig('./results/{}/clf_ALL.png'.format(dataset))



print("Best for COTAN: {}".format(best_for_COTAN))
print("Worst for COTAN: {}".format(worst_for_COTAN))
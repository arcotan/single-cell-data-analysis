import pandas as pd

# https://www.tutorialspoint.com/matplotlib/matplotlib_bar_plot.htm


"""DATASET_TAGS = {
	'tabula-muris-heart': 'tabula muris heart',
	'tabula-muris-marrow_P7_3': 'tabula muris marrow',
	'peripheal-blood': 'peripheal blood',
	'zheng-4': 'zheng 4',
	'zheng-8': 'zheng 8'
}"""

DATASET_TAGS = {
	'tabula-muris-heart': 'tabula muris heart',
	#'tabula-muris-marrow_P7_3': 'tabula muris marrow',
	'peripheal-blood': 'peripheal blood'#,
	#'zheng-4': 'zheng 4',
	#'zheng-8': 'zheng 8'
}

DATASET_ORDER = {
	'tabula-muris-heart': 1,
	'tabula-muris-marrow_P7_3': 2,
	'peripheal-blood': 3,
	'zheng-4': 4,
	'zheng-8': 5
}

TOOL_TAGS = {
    'monocle':'Monocle',
    'scanpy': 'Scanpy',
    'seurat': 'Seurat',
    'scvi': 'scvi-tools',
    'COTAN': 'COTAN'
}

path = './results/aggregate/{}/scores.csv'

scores_df = pd.DataFrame(columns=['accuracy', 'entropy', 'purity', 'silhouette', 'tool', 'dataset'])
for dataset in DATASET_TAGS.keys():
	file_path = path.format(dataset)
	df = pd.read_csv(file_path)
	df['dataset'] = dataset
	scores_df = pd.concat([scores_df, df])
	
scores_df['dataset_pos'] = scores_df['dataset']
scores_df = scores_df.replace({'dataset_pos' : DATASET_ORDER})

##### plot

import matplotlib.pyplot as plt
import numpy as np

SCORE = 'entropy'
datasets = ('tabula muris heart', 'tabula muris marrow', 'peripheal blood', 'zheng 4', 'zheng 8')
#penguin_means = {
#    'cotan': (18.35, 18.43, 14.98), # list of 5 scores (one for each dataset in order)
#    'seurat': (38.79, 48.83, 47.50),
#    'monocle' : (189.95, 195.82, 217.19),
#    'scvitool' : (189.95, 195.82, 217.19),
#    'scanpy' : (189.95, 195.82, 217.19),
#}
penguin_means = {}
for tool in TOOL_TAGS.keys():
	tool_scores = scores_df[scores_df['tool'] == tool]
	tool_scores.sort_values(by='dataset_pos')
	penguin_means[tool] = (score for score in tool_scores[SCORE])

print(penguin_means)

x = np.arange(len(datasets))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in penguin_means.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Length (mm)')
ax.set_title('Penguin attributes by species')
ax.set_xticks(x + width, datasets)
ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 250)

plt.show()
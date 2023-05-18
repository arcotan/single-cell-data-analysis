import pandas as pd
import plotly.graph_objects as go
import os

DATA_DIR = './results/aggregate'

# DATASET_TAGS = ['tabula-muris-heart', 'tabula-muris-marrow_P7_3', 'peripheal-blood', 'kumar-4-hard', 'kumar-8-hard']
DATASET_TAGS = ['tabula-muris-heart']

colors = [
    '#1f77b4',  # muted blue
    '#ff7f0e',  # safety orange
    '#2ca02c',  # cooked asparagus green
    '#d62728',  # brick red
    '#9467bd',  # muted purple
    '#8c564b',  # chestnut brown
    '#e377c2',  # raspberry yogurt pink
    '#7f7f7f',  # middle gray
    '#bcbd22',  # curry yellow-green
    '#17becf'   # blue-teal
]

def plot_sankey(labels, source, target, value, title):
    fig = go.Figure(data=[go.Sankey(
        node = dict(
        pad = 15,
        thickness = 20,
        line = dict(color = "black", width = 0.5),
        label = labels,
        color = colors[0:len(source)-1] + colors[0:len(target)-1]
        ),
        link = dict(
        source = source,
        target = target,
        value = value
        ))])

    fig.update_layout(title_text=title, font_size=10)
    fig.show()

for dataset in DATASET_TAGS:
    cur_path = DATA_DIR + '/' + dataset + '/'
    # check if labels.csv exists
    if os.path.exists(cur_path + "labels.csv"):
        labels = pd.read_csv(cur_path + "labels.csv")

        if labels.columns[-1] != 'true_labels':
            print("ERROR: true_label not found for dataset " + dataset)
        else:
            tool_tags = labels.columns[2:-1]
            cluster_num = labels['true_labels'].nunique()
            plot_label = []
            for i in range(cluster_num):
                plot_label.append('cluster' + str(i))
            for i in range(cluster_num):
                mapping_path = "./dataset/" + dataset + "-filtered/mapping.csv"
                mapping_df = pd.read_csv(mapping_path)
                go_id = mapping_df[mapping_df['id'] == i]['go'].values
                plot_label.append(go_id)

            for tool in tool_tags:
                # generate confusion matrix between labels
                confusion_matrix = pd.crosstab(labels[tool], labels['true_labels'], colnames=['Predicted'], rownames=['True'], margins=True)
                cur_source = []
                cur_target = []
                cur_value = []
                for i in range(cluster_num):
                    for j in range(cluster_num):
                        if confusion_matrix.iloc[i, j] != 0:
                            cur_source.append(i)
                            cur_target.append(cluster_num + j)
                            cur_value.append(confusion_matrix.iloc[i, j])
                print(cur_source, cur_target, cur_value)
                plot_sankey(plot_label, cur_source, cur_target, cur_value, tool + ' ' + dataset)
            

"""fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = ["A1", "A2", "B1", "B2", "C1", "C2"],
      color = "blue"
    ),
    link = dict(
      source = [0, 1, 0, 2, 3, 3], # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = [2, 3, 3, 4, 4, 5],
      value = [8, 4, 2, 8, 4, 2]
  ))])

fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
fig.show()"""
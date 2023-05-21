echo 'Starting Clustering Analysis'
chmod +x coderun.sh
./coderun.sh
R CMD BATCH ./results/Clustering_Analysis.R
echo 'End Clustering Analysis'

echo 'Starting Num of Cells Genes Check'
ipython -c "%run num_cells_genes.ipynb"
echo 'End Num of Cells Genes Check'

echo 'Starting Classifier'
python3 clf.py
echo 'End Classifier'

echo 'Starting Classifier Plots'
ipython -c "%run clf_plots.ipynb"
echo 'End Classifier Plots'

echo 'Starting Classifier Histograms'
ipython -c "%run classifier_histograms.ipynb"
echo 'End Classifier Histograms'

echo 'Starting Sankey Plots'
ipython -c "%run sankey.ipynb"
echo 'End Sankey Plots'

echo 'Starting tabula muris heart markers'
python3 tabula-muris-heart-markers.py
end 'End tabula muris heart markers'

echo 'Starting tabula muris heart markers'
ipython -c "%run tabula-muris-heart-markers-hist.ipynb"
end 'End tabula muris heart markers'

echo 'Starting Cluster Scores Plots'
ipython -c "%run cluster_scores_plot.ipynb"
echo 'End Cluster Scores Plots'

exit

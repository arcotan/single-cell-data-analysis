echo 'Starting Clustering Analysis'
chmod +x coderun.sh
./coderun.sh
R CMD BATCH ./analysis.R
echo 'End Clustering Analysis'

echo 'Starting Num of Cells Genes Check'
ipython -c "%run ./scripts/num_cells_genes.ipynb"
echo 'End Num of Cells Genes Check'

echo 'Starting Classifier'
python3 ./scripts/clf.py
echo 'End Classifier'

echo 'Starting Classifier Plots'
ipython -c "%run ./plots/clf_plots.ipynb"
echo 'End Classifier Plots'

echo 'Starting Classifier Histograms'
ipython -c "%run ./plots/classifier_histograms.ipynb"
echo 'End Classifier Histograms'

echo 'Starting Sankey Plots'
ipython -c "%run ./plots/sankey.ipynb"
echo 'End Sankey Plots'

echo 'Starting tabula muris heart markers'
python3 ./scripts/tabula-muris-heart-markers.py
end 'End tabula muris heart markers'

echo 'Starting tabula muris heart markers'
ipython -c "%run ./plots/tabula-muris-heart-markers-hist.ipynb"
end 'End tabula muris heart markers'

echo 'Starting Cluster Scores Plots'
ipython -c "%run ./plots/cluster_scores_plot.ipynb"
echo 'End Cluster Scores Plots'

exit

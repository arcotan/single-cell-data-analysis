# TO RUN THIS SCRIPT ALL THE TOOLS MUST BE ALREADY RUNNED ON EACH DATASET
# This script runs all the analysis and plots for the project
echo 'Starting Clustering Analysis'
((sleep 2 && code ./analysis.Rout) & 
(R CMD BATCH ./scripts/analysis.R))
echo 'End Clustering Analysis'

echo 'Starting Num of Cells Genes Check'
cd ./scripts
ipython -c "%run ./dataset_sizes.ipynb"
echo 'End Num of Cells Genes Check'

echo 'Starting Classifier'
python3 ./clf.py
cd ..
echo 'End Classifier'

echo 'Starting Classifier Plots'
cd ./plots
ipython -c "%run ./clf_curves.ipynb"
echo 'End Classifier Plots'

echo 'Starting Classifier Histograms'
ipython -c "%run ./clf_histogram.ipynb False"
echo 'End Classifier Histograms' 

echo 'Starting Sankey Plots'
ipython -c "%run ./sankey.ipynb False"
cd ..
echo 'End Sankey Plots'

echo 'Starting tabula muris heart markers script'
cd ./scripts
python3 ./tabula-muris-heart-markers.py
cd ..
echo 'End tabula muris heart markers'

echo 'Starting tabula muris heart markers notebook'
cd ./plots
ipython -c "%run ./tabula-muris-heart-markers-hist.ipynb False"
echo 'End tabula muris heart markers'

echo 'Starting Cluster Scores Plots'
ipython -c "%run ./cluster_scores_plot.ipynb False"
echo 'End Cluster Scores Plots'

echo 'Starting Enrichment Scores Hist'
ipython -c "%run ./enrichment_scores_hist.ipynb False"
cd ..
echo 'End Enrichment Scores Hist'

exit

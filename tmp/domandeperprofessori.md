# tante belle domande su cotan e molto altro

- Tabula muris ha geni mitocondriali? Cercando "Mthfd1l" su ontology appare qualcosa di simile ma
qui dicono che non sono stati sequenziali https://github.com/czbiohub/tabula-muris/issues/221
- 

# COSE

- quale database usare per i geni di peripheral blood
- monocle può usare leiden (con risoluzione, usat anche da altri tool) o louvian (senza risoluzione)
- parametri usati dagli autori di tabula muris per analisi https://github.com/czbiohub/tabula-muris/tree/master/00_data_ingest/02_tissue_analysis_rmd
- seurat usa AUROC: https://satijalab.org/seurat/reference/findmarkers
- Invece di enrichment a mano monocle offre Garnett: https://cole-trapnell-lab.github.io/monocle3/docs/clustering/ e Seurat https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html#cell-type-annotation-using-singler, possiamo valutare # cellule correttamente identificate su tabula muris

# PLOT

- tipo chord diagram per casi particolari (altrimenti sono 25) / confusion matrix
- venn diagram (intersezioni top genes per ogni cluster / intersezioni top 120 genes)
- confronto clustering
- box plot con espressione (controllare tanto espressi)
- heatmap con markers
- dot plot https://satijalab.org/seurat/articles/visualization_vignette.html 
   https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html
- violin plot 1 vs rest per vedere distribuzioni dei markers nei clusters


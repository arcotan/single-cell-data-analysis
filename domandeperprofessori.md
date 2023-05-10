# tante belle domande su cotan e molto altro

- come fare enrichment (come usare david, specie di mouse)
- il preprocessing deve essere diverso con cotan? Ad esempio ha senso togliere i geni poco espressi?
- significato dei plot
- perche' prima de e poi clustering?
- correlazioni tra gene.id, cell ontology class e subsetA (riferimento al file analisi), per allenare il modello eventualmente
- mostrare le pipeline
- usare gli stessi parametri o quelli standard? altrimenti diventano uno la copia dell'altro (esempio scanpy e seurat)
- quale test statistico utilizzare?

- quali topi prendere? (10X_P7_2 e/o 10X_P7_3)
- come calcolare silhouette (nel latent space di PCA? articolo paper: https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue    codice: https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36)
- controllare se bechmark 3 ha labels (se non ce l'ha non sappiamo nemmeno k)
- top 20 è okay? +/-?
- se un gene è un marcatore per un cluster ma è mediamente abbastanza espresso anche negli altri non ha significato bilogico? viceversa se poco espresso?
- DE con funzioni particolari di tool diversi dati true labels o DE con stesso test dati risultati clustering?

# COSE

- parametri usati dagli autori di tabula muris per analisi https://github.com/czbiohub/tabula-muris/tree/master/00_data_ingest/02_tissue_analysis_rmd

# PLOT

- tipo chord diagram per casi particolari (altrimenti sono 25) / confusion matrix
- venn diagram (intersezioni top genes per ogni cluster / intersezioni top 120 genes)
- confronto clustering
- box plot con espressione (controllare tanto espressi)
- heatmap con markers
- dot plot https://satijalab.org/seurat/articles/visualization_vignette.html 
   https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html
- violin plot 1 vs rest per vedere distribuzioni dei markers nei clusters


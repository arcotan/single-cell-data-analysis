# cose da fare

- cercare algoritmi specifici per ogni tool (resta da vedere monocle)
- classificazione (implementare misura del ranking)
- come organizzare il report?
- enrichment analysis e venn diagrams
- valutare la presenza dei marcatori che hanno usato in tabula muris per etichettare cellule nel ranking prodotto dalla DE
- fare plot De per tool in python
- rimuovere il filtering dallo script di seurat 
- come calcolare silhouette (nel latent space di PCA? articolo paper: https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue    codice: https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36)

- seurat usa PCA con 10 variables features, gli altri?


- uniformare letture dataset in modo che tutti gli script leggano unica matrice di espressione, lo script che crea questa matrice deve essere specifico per ogni dataset. Eventualmente (ma forse non è una buona idea) inclucere nella nuova matrice di espressione anche degli step che adesso si fanno separatamente come ad esempio la normalizzazione ed il logaritmo.
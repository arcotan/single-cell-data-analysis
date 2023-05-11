# cose da fare

- capire se in peripheal blood le 19 labels individuate hanno senso (il processo che abbiamo usato è okay)
- classificazione (implementare misura del ranking, calcolo del numero di marcatori da usare)
- come organizzare il report?
- ERCC
- perchè genes.tsv ci sono due nomi
- enrichment analysis
- capire se scalare su scanpy e scivitools (meglio con)
- valutare la presenza dei marcatori che hanno usato in tabula muris per etichettare cellule nel ranking prodotto dalla DE
- fare plot De per tool in python
- rimuovere il filtering dallo script di seurat 
- come calcolare silhouette (nel latent space di PCA? articolo paper: https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue    codice: https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36)

- seurat usa PCA con 10 variables features, gli altri?


- Eventualmente (ma forse non è una buona idea) inclucere nella nuova matrice di espressione anche degli step che adesso si fanno separatamente come ad esempio la normalizzazione ed il logaritmo.
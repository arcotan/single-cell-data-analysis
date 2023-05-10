# cose da fare

- cercare algoritmi specifici per ogni tool
- classificazione
- come organizzare il report?
- enrichment analysis?
- uniformare letture dataset in modo che tutti gli script leggano unica matrice di espressione, lo script che crea questa matrice deve essere specifico per ogni dataset. Eventualmente (ma forse non è una buona idea) inclucere nella nuova matrice di espressione anche degli step che adesso si fanno separatamente come ad esempio la normalizzazione ed il logaritmo.
- valutare la presenza dei marcatori che hanno usato in tabula muris per etichettare cellule nel ranking prodotto dalla DE
- fare plot De per tool in python
- decidere come valutare ranking nella cassificazione? anche in tabula muris fanno qualcosa di simile
- seurat usa AUROC: https://satijalab.org/seurat/reference/findmarkers
- Invece di enrichment a mano monocle offre Garnett: https://cole-trapnell-lab.github.io/monocle3/docs/clustering/ e Seurat https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html#cell-type-annotation-using-singler, possiamo valutare # cellule correttamente identificate su tabula muris
- cambiare png in eps o pdf


# future

- applicare ad altri dataset

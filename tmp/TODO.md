Ordine e labels: COTAN, Seurat, Scanpy, Monocle, scvi-tools
Ordine dataset: Tabula Muris Heart, Marrow, Peripheral, Zheng 4 e 8

PLOT CHE ABBIAMO:
- istogrammi (scores clustering, classificazione f1, arricchimento, tabula muris check)
- venn
- sankey (stessi colori pca)
- enrichr (va bene colore usato)
- pca (arcobaleno)
- classificatore curve con feature crescenti

TODO:
- titoli ai plot (nomi dei dataset)
- nomi legende: dataset tag e tools? (venn)
- colori del Set2
- Scorriamo notebooks per filtering e pipeline specifiche
- Citiamo cosa abbiamo fatto per ricavare le labels di peripheral
- Plot PCA o altro dim reduction?
- Confrontare plot sito tabula muris?
- Slide con riferimenti ai tool alla fine
- Allineamento clustering (limitazioni e assunzioni), definizioni metriche
- 5 slides con: un plot PCA del ground truth + 5 sankey per ogni tool, poi decidiamo poi quale tenere
- Istogrammi con scores
- Citare plot con DE dei top genes (scorriamoli e vediamo se ci sono cose strane)
- Diagrammi di venn e enrichr di esempio!
- Diagrammi di venn in cui colori solo centro / solo fuori / tutto tranne dentro con colori diversi
- **Istogramma con quanti cluster sono giusti con arricchimento**
- **Tabula Muris enrichr validazione markers, istogramma con quanti marcatori veri compaiono nella top 50/40/30…**
- 5 plot della classificazione con le curve
- Istogramma della classificazione
- NA cotan?
- Tempi computazione
- Punti forza e debolezza di ciascun tool
- Capire meglio se la tabella sui tools è giusta

- Slides aggiuntive da mettere alla fine se capitano domande

TABELLA ARRICCHIMENTO

Enrichment di markers di ciascun tool - intersezione comune a tutti
Dataset, #cluster, seurat, monocle, COTAN, 
tabula-muris-heart, 5, 4, 3, 5

Di quanti cluster l'arricchimento corrispondeva alla gene ontology?

TABULA MURIS HEART MARKER VALIDATION

Quanti marcatori del cluster compaiono nella top 50

Cluster,#marcatori,Seurat,Monocle,...
fibroblast,3,?,?
endothelial,
...
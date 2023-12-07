# Number of clusters found with default parameters
| dataset| COTAN | scanpy | scvi-tools | seurat | monocle | celltypist |
| :----: | :---: | :----: | :--------: | :----: | :-----: | :--------: |
| PBMC1 | 28 | 18 |  12 | 11 | 3 | 15 |
| PBMC2 | 31 | 19 | 16 | 16 | 3 | 16 |
| PBMC3 | 61 | 22 | 19 | 18 | 4 | 19 |
| PBMC4 | 43 | 22 | 16 | 19 | 3 | 17 |

# Labeling with Antibody derived tags

## PBMC1
![image](./dataset/PBMC1-Filtered/10x/heatmap.png)
![image](./plots/pbmc1_sankey.png)

| celltypist cluster id | celltypist cluster identity |
| :-------------------: | :-------------------------- |
| 0 | Classical monocytes |
| 1 | Tcm/Naive helper T cells |
| 2 | Regulatory T cells |
| 3 | DC2 |
| 4 | Tem/Effector helper T cells |
| 5 | Naive B cells |
| 6 | MAIT cells |
| 7 | CD16+ NK cells |
| 8 | Tcm/Naive cytotoxic T cells |
| 9 | Non-classical monocytes |
| 10 | Tem/Temra cytotoxic T cells |
| 11 | Tem/Trm cytotoxic T cells |
| 12 | pDC |
| 13 | Memory B cells |
| 14 | Megakaryocytes/platelets |

### Scores
| accuracy | entropy | purity | silhouette | NMI | Adj_Rand_Index | tool |
| :------: | :-----: | :----: | :--------: | :-: | :------------: | :--: |
| 0.534072022160665 | 0.468433704758432 | 0.577285318559557 | 0.412082723163031 | 0.576502096293195 | 0.953268085545066 | monocle |
| 0.674318507890961 | 0.176770424028319 | 0.83751793400287 | 0.134876401838029 | 0.763087788644786 | 0.559873236269802 | scanpy |
| 0.390027700831025 | 0.171131878500729 | 0.838504155124654 | 0.183325293692934 | 0.794971537394654 | 0.67103279528902 | seurat |
| 0.270360110803324 | 0.171754199850388 | 0.804986149584488 | 0.0838710472015484 | 0.75245330509039 | 0.583024251768234 | scvi-tools |
| 0.75448192038894 | 0.134446726354925 | 0.878152537222729 | 0.148614434893049 | 0.813650352649089 | 0.730101427112584 | COTAN |

## PBMC2
![image](./dataset/PBMC2-Filtered/10x/heatmap.png)
![image](./plots/pbmc2_sankey.png)

| celltypist cluster id | celltypist cluster identity |
| :-------------------: | :-------------------------- |
| 1 | Tem/Trm cytotoxic T cells |
| 2 | Tem/Effector helper T cells |
| 3 | Tcm/Naive helper T cells |
| 4 | CD16+ NK cells |
| 5 | MAIT cells |
| 6 | DC2 |
| 7 | Classical monocytes |
| 8 | Tcm/Naive cytotoxic T cells |
| 9 | Non-classical monocytes |
| 10 | CD16- NK cells |
| 11 | Naive B cells |
| 12 | Memory B cells |
| 13 | pDC |
| 14 | Regulatory T cells |
| 15 | HSC/MPP |
| 16 | Megakaryocytes/platelets |

### Scores
| accuracy | entropy | purity | silhouette | NMI | Adj_Rand_Index | tool |
| :------: | :-----: | :----: | :--------: | :-: | :------------: | :--: |
| 0.0328958298720718 | 0.511893876590759 | 0.513042033560392 | 0.447336360694943 | 0.48310954481693 | 0.748608157689483 | monocle | 
| 0.670985915492958 | 0.173097979852628 | 0.797183098591549 | 0.135114211594495 | 0.729596971527379 | 0.509794938339751 | scanpy |
| 0.734839674364512 | 0.139250356205862 | 0.847482970593122 | 0.160773352731945 | 0.755414119223764 | 0.522141016645438 | seurat |
| 0.68981558398405 | 0.16539230446373 | 0.824223292905798 | 0.137456937702708 | 0.721636369909799 | 0.48099500276417 | scvi-tools |
 0.822663395047959 | 0.106734440894982 | 0.87798349319652 | 0.216846167658323 | 0.819508774894713 | 0.658310337129159 | COTAN |

## PBMC3
![image](./dataset/PBMC3-Filtered/10x/heatmap.png)
![image](./plots/pbmc3_sankey.png)

| celltypist cluster id | celltypist cluster identity |
| :-------------------: | :-------------------------- |
| 1 | Tcm/Naive helper T cells |
| 2 | Classical monocytes |
| 3 | Memory B cells |
| 4 | Tcm/Naive cytotoxic T cells |
| 5 | CD16+ NK cells |
| 6 | Tem/Effector helper T cells |
| 7 | DC2 |
| 8 | MAIT cells |
| 9 | Naive B cells |
| 10 | Non-classical monocytes |
| 11 | Tem/Trm cytotoxic T cells |
| 12 | Plasma cells |
| 13 | Regulatory T cells |
| 14 | NK cells 
| 15 | Megakaryocytes/platelets |
| 16 | Epithelial cells |
| 17 | DC1 |
| 18 | pDC |
| 19 | HSC/MPP |

### Scores
| accuracy | entropy | purity | silhouette | NMI | Adj_Rand_Index | tool |
| :------: | :-----: | :----: | :--------: | :-: | :------------: | :--: |
| 0.475235849056604 | 0.423986386328212 | 0.569212626995646 | 0.362639216723731 | 0.60083021907912 | 0.762309381636632 | monocle |
| 0.663942401245379 | 0.161577294576172 | 0.817960692741779 | 0.126767711592436 | 0.736294077113439 | 0.556489492175471 | scanpy | 
| 0.751088534107402 | 0.122291676762917 | 0.873185776487663 | 0.146239000944664 | 0.783918806534948 | 0.606468875865081 | seurat |
| 0.717343976777939 | 0.173679090632666 | 0.79617198838897 | 0.147064633022423 | 0.728299366408049 | 0.577217688322832 | scvi-tools |
| 0.943786982248521 | 0.0489899209299561 | 0.96560650887574 | 0.260378965397053 | 0.921155209292118 | 0.92592979007967 | COTAN |


## PBMC4
![image](./dataset/PBMC4-Filtered/10x/heatmap.png)
![image](./plots/pbmc4_sankey.png)

| celltypist cluster id | celltypist cluster identity |
| :-------------------: | :-------------------------- |
| 1 | MAIT cells
| 2 | Naive B cells
| 3 | Tcm/Naive helper T cells
| 4 | Tem/Temra cytotoxic T cells
| 5 | Classical monocytes
| 6 | DC2
| 7 | Tcm/Naive cytotoxic T cells
| 8 | Tem/Effector helper T cells
| 9 | CD16+ NK cells
| 10 | Non-classical monocytes
| 11 | Memory B cells
| 12 | Tem/Trm cytotoxic T cells
| 13 | Megakaryocytes/platelets
| 14 | Regulatory T cells
| 15 | pDC
| 16 | Plasma cells
| 17 | ILC3

### Scores
| accuracy | entropy | purity | silhouette | NMI | Adj_Rand_Index | tool |
| :------: | :-----: | :----: | :--------: | :-: | :------------: | :--: | 
| 0.000980804259492784 | 0.412554830187043 | 0.620568866470506 | 0.344929310461403 | 0.616374371955131 | 0.807514770527812 | monocle |
| 0.680901542111507 | 0.13928737712969 | 0.848669716997119 | 0.14373826389226 | 0.771759026594095 | 0.578919064651214 | scanpy |
| 0.66297558205565 | 0.114281272994342 | 0.867830777967064 | 0.156123561272193 | 0.768775901756073 | 0.502952478483863 | seurat |
| 0.627854840969595 | 0.15482862855921 | 0.837326607818411 | 0.147943469695352 | 0.748858810993074 | 0.51797069778846 | scvi-tools |
| 0.82691131498471 | 0.12357738164552 | 0.857900101936799 | 0.23366038135263 | 0.830215446325249 | 0.704023119079152 | COTAN |
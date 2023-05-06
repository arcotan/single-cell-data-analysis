library(COTAN)
library(zeallot)
library(ggplot2)
library(factoextra)
library(qpdf)


source("utils.R")

outDir <- tempdir()

IN_DATA_DIR = "./dataset/tabulamuris/droplet/Heart_and_Aorta-10X_P7_4"
IN_LABEL_DIR = "./dataset/tabulamuris"
CHANNEL = "10X_P7_4"

experiment_data = load_data(IN_DATA_DIR, IN_LABEL_DIR, CHANNEL)
data <- experiment_data$data

obj = COTAN(raw = data)
obj = initializeMetaDataset(obj,
                            GEO = "-",
                            sequencingMethod = "10x",
                            sampleCondition = "Heart_and_Aorta")

ECDPlot(obj, yCut = 700) # a cosa serve
cellSizePlot(obj) # per ogni cellula la somma dei counts

mit <- mitochondrialPercentagePlot(obj, genePrefix = "^Mt")
mit[["plot"]]

cells_to_rem <- getCells(obj)[getCellsSize(obj) > 15000] # la somma dei counts di ciascuna cellula deve essere < 15000 (elimina doublets)
obj <- dropGenesCells(obj, cells = cells_to_rem)

cellGeneNumber <- sort(colSums(as.data.frame(getRawData(obj) > 0)), decreasing = FALSE)
cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > 4000] # outliers
obj <- dropGenesCells(obj, cells = cells_to_rem)
genesSizePlot(obj) # per ogni cellula il numero di reads > 0

to_rem <- mit[["sizes"]][["mit.percentage"]] > 0.8
cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
obj <- dropGenesCells(obj, cells = cells_to_rem)

mit <- mitochondrialPercentagePlot(obj, genePrefix = "^Mt")
mit[["plot"]]

obj <- clean(obj)
c(pcaCellsPlot, pcaCellsData, genesPlot, UDEPlot, nuPlot) %<-% cleanPlots(obj)
pcaCellsPlot # ?
genesPlot
UDEPlot
nuPlot

ggplot(subset(pcaCellsData[pcaCellsData$PC1 < 0,], groups == "A"),
                       aes(x = PC1, y = PC2, colour = groups)) +
  geom_point(alpha = 0.5, size = 3L) +
  geom_point(data = subset(pcaCellsData[pcaCellsData$PC1 < 0,], groups != "A"),
             aes(x = PC1, y = PC2, colour = groups),
             alpha = 0.8, size = 3L) +
  scale_color_manual(groups, values = c("A" = "#8491B4B2",
                                        "B" = "#E64B35FF")) +
  plotTheme(
    "pca")

obj <- clean(obj)
cells_to_rem <- rownames(pcaCellsData)[pcaCellsData[["groups"]] == "B"]
obj <- dropGenesCells(obj, cells = cells_to_rem)

UDEPlot
nuPlot

nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
yset = 0.20 # the threshold to remove low UDE cells
plot.ude <- ggplot(nuDf, aes(x = n, y = nu)) +
  geom_point(colour = "#8491B4B2", size = 1) +
  xlim(0, 400) +
  ylim(0,   1) +
  geom_hline(yintercept = yset, linetype = "dashed",
             color = "darkred") +
  annotate(geom = "text", x = 200, y = 0.25, 
           label = paste0("to remove cells with nu < ", yset), 
           color = "darkred", size = 4.5)

plot.ude

obj <- addElementToMetaDataset(obj, "Threshold low UDE cells:", yset) 

cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < yset]
obj <- dropGenesCells(obj, cells = cells_to_rem)


obj <- clean(obj)
c(pcaCellsPlot, pcaCellsData, genesPlot, UDEPlot, nuPlot) %<-% cleanPlots(obj)

pcaCellsPlot

#EVVIVA
obj = estimateDispersionBisection(obj, cores = 10)
obj <- calculateCoex(obj)
saveRDS(obj, file = file.path("./dataset/tabulamuris/", paste0("Heart_and_Aorta", ".cotan.RDS")))

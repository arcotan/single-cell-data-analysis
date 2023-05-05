library(COTAN)

DATA_PATH = "./filtered_dataset/tabulamuris/data_10X_P7_4.csv"

outDir <- tempdir()

data <- read.csv(DATA_PATH, row.names = 1L)

# bugged
# obj <- automaticCOTANObjectCreation(
#   raw = data,
#   GEO = "-",
#   sequencingMethod = "10x",
#   sampleCondition = "Heart_and_Aorta",
#   saveObj = TRUE,
#   outDir = outDir,
# )

obj = COTAN(raw = data)
obj = initializeMetaDataset(obj,
                            GEO = "-",
                            sequencingMethod = "10x",
                            sampleCondition = "Heart_and_Aorta")

logThis(paste("n cells", getNumCells(obj)), logLevel = 1)

# efficiency filtering skipped (todo?)

obj = estimateDispersionBisection(obj)
obj <- calculateCoex(obj)

quant.p = calculateGDI(obj)
head(quant.p)
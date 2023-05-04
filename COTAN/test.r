library(COTAN)

DATA_PATH = "./filtered_dataset/tabulamuris/data_10X_P7_4.csv"

outDir <- tempdir()

data <- read.csv(DATA_PATH, row.names = 1L)

obj <- automaticCOTANObjectCreation( # drops genes!
  raw = data,
  GEO = "-",
  sequencingMethod = "10x",
  sampleCondition = "Heart_and_Aorta",
  saveObj = TRUE,
  outDir = outDir,
)
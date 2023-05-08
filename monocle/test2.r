rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(reticulate))
suppressMessages(library(devtools))
suppressMessages(library(monocle))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))


DATA_DIR = "./dataset/tabulamuris/droplet/Heart_and_Aorta-10X_P7_4"
LABEL_DIR = "./dataset/tabulamuris"
CHANNEL = "10X_P7_4"


pulm <- read.csv(DATA_DIR)

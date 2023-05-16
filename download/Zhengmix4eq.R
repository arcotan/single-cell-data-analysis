library(ExperimentHub)
library(SingleCellExperiment)

DATASET_DIR = "./dataset/zheng-4"

eh <- ExperimentHub()
q <- query(eh, "DuoClustering2018")
q$title[grep("zheng", q$title, ignore.case=TRUE)]
id_zheng4 = q$ah_id[q$title == "sce_full_Zhengmix4eq"]

db = q[[id_zheng4]]

if (!dir.exists(DATASET_DIR)) {
	dir.create(DATASET_DIR)
}

saveRDS(db, file.path(DATASET_DIR, "zheng-4.rds"))



library(ExperimentHub)
library(SingleCellExperiment)

DATASET_DIR = "./dataset/zheng-8"

eh <- ExperimentHub()
q <- query(eh, "DuoClustering2018")
q$title[grep("zheng", q$title, ignore.case=TRUE)]
id_zheng8 = q$ah_id[q$title == "sce_full_Zhengmix8eq"]

db = q[[id_zheng8]]

if (!dir.exists(DATASET_DIR)) {
	dir.create(DATASET_DIR)
}

saveRDS(db, file.path(DATASET_DIR, "sce_full_Zhengmix8eq.rds"))



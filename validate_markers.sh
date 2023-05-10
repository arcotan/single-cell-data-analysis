#!/bin/bash

# ---- TABULA MURIS ----

# Seurat
# TODO: change paths
python ./markers_validation.py \
--data-path "./filtered_dataset/tabulamuris/data_10X_P7_4/" \
--markers-path "./results/markers_10X_P7_4_seurat.csv" \
--labels-path "./filtered_dataset/tabulamuris/labels_10X_P7_4.csv" \
--out-path "./results/" \

# Scanpy
# ...

# ---- KUMAR 4 HARD ----
# ...
# --simulated
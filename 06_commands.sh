#Work in Progress
# Written by Christos Vlachos, Markus Boesch and Emiel Geeraerts
Rscript 01_pre_process.R --out-folder="results" --raw-matrix="raw_matrices" --matrix="matrices"
Rscript 02_rename_cells.R --out-folder="results" --rds="results/objects/preprocessed_list.rds.gz"
Rscript 03_QC.R --out-folder="results" --rds="results/objects/renamed_objects.rds.gz" --QC="300,10000,0,999999,50,99999" #set QC based min_feauture,max_feature,min_count, max_count, MT%, RB%
Rscript 04_dimension_reduction.R --out-folder="results" --rds="results/objects/QC_list.rds.gz" --norm-method="SCT" --regress-out="percent.mt" --exclude="merged" #exclude the merged sample from integration and others if needed
Rscript 05_integration.R --out-folder="results" --rds="results/objects/Normalized_merged_curated_final.rds.gz"



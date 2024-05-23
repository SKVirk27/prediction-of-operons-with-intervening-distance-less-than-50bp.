Operon Prediction in Bacterial Genomes
Description: This script predicts operons in bacterial genomes by calculating intergenic distances between genes and identifying pairs of genes likely to form operons.

Key Features:

Calculates intergenic distances between gene pairs.
Identifies operons based on defined criteria (e.g., intergenic distance < 50bp).
Outputs include B_subtilis_OperonPrediction.csv, E_coli_OperonPrediction.csv, Halob_NRC1_OperonPrediction.csv, Synec_PCC_OperonPrediction.csv, Metagenome_OperonPrediction.csv.
Required Libraries: os, pandas, gffpandas

Execution:

Ensure the required files (B_subtilis_168.ptt.gz, E_coli_K12_MG1655.ptt.gz, Halobacterium_NRC1.ptt.gz, Synechocystis_PCC6803_uid159873.ptt.gz, 2088090036 (1).gff, main.py) are in the same folder.
Open the project in PyCharm and set the working folder as Assignment3.
Run main.py to execute the script and predict operons in the bacterial genomes.

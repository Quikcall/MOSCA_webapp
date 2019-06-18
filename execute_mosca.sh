#!/bin/sh
python MOSCA/scripts/mosca.py	--files	/media/sf_shared_folder/MOSCA_app/templates/layout.html:/media/sf_shared_folder/MOSCA_app/templates/login.html	--sequencing-technology paired	--assembler	metaspades	--annotation-database	/mosca/Databases/annotation_databases/uniprot.fasta	--output	/home/mario/shared_folder/MOSCA_app/MOSCA/Output	--output-level	min	--type-of-data	metatranscriptomics 	--conditions Etanol	--threads	4	--memory	10	--marker-gene-set	40 
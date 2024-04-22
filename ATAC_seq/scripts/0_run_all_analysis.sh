# Process all the analysis for the ATAC-seq data

# preprocess the data for the analysis
bash ATAC_seq/scripts/1_preprocess.sh

# run the differential accessibility analysis and generate the plots
CONTRASTS="hsc_cfue hsc_ery hsc_cmp cmp_ery cfue_ery cmp_cfue"
for contrast in $CONTRASTS
do
	Rscript ATAC_seq/scripts/2_differential_accessibility_analysis.R $contrast
	python ATAC_seq/scripts/3_chromatin_accessibility_vs_gene_expression.py --contrast $contrast
	python ATAC_seq/scripts/4_plot_GO_enrichment.py --contrast $contrast
done

# perform clustering analysis
Rscript ATAC_seq/scripts/5_clustering_analysis.R
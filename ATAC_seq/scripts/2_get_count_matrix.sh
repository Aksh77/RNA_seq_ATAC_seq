# make a list of pairs
CONTRASTS="hsc_cfue hsc_erythroblast hsc_cmp cmp_erythroblast cfue_erythroblast cmp_cfue"
for contrast in $CONTRASTS
do
    # get cell lines for each contrast
    cell1=$(echo $contrast | cut -d'_' -f1)
	cell2=$(echo $contrast | cut -d'_' -f2)	
	echo $cell1 $cell2

	# merge consensus peaks into peakset
	peaks1=ATAC_seq/data/peaks/${cell1}/${cell1}_peaks.narrowPeak
	peaks2=ATAC_seq/data/peaks/${cell2}/${cell2}_peaks.narrowPeak
	merged_peakset=ATAC_seq/data/peaks/${cell1}_${cell2}_peaks.narrowPeak
	bedops -u $peaks1 $peaks2 > $merged_peakset

	# get files for featureCounts
	output_dir=ATAC_seq/data/output_data/featureCounts/${contrast}
	mkdir -p $output_dir
	contrast_bed=${output_dir}/${contrast}_peaks.bed
	contrast_saf=${output_dir}/${contrast}_peaks.saf
	contrast_counts=${output_dir}/${contrast}_peaks.counts
	contrast_counts_tsv=${output_dir}/${contrast}_peaks.tsv
	awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="'${contrast}'_merged_"++nr;  print peakid,$1,$2,$3,"."}' $merged_peakset > $contrast_saf
	awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="'${contrast}'_merged_"++nr;  print $1,$2,$3,peakid,"0","."}' $merged_peakset > $contrast_bed

	# get BAM files for the cell lines
	BAM_DIR=ATAC_seq/data/input_data
	BAM1=$(ls $BAM_DIR/${cell1}/*.bam | head -n 1)
	BAM2=$(ls $BAM_DIR/${cell1}/*.bam | tail -n 1)
	BAM3=$(ls $BAM_DIR/${cell2}/*.bam | head -n 1)
	BAM4=$(ls $BAM_DIR/${cell2}/*.bam | tail -n 1)
	echo $BAM1 $BAM2 $BAM3 $BAM4

	# get the count matrix for peaks
	featureCounts -F SAF -a $contrast_saf --fracOverlap 0.2 -o $contrast_counts $BAM1 $BAM2 $BAM3 $BAM4

	# remove first line starting with #
	awk '(NR>1)' $contrast_counts > $contrast_counts_tsv

done

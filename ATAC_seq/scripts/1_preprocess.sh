# preprocess the data for the analysis
for bigbed in $(ls ATAC_seq/data/input_data/*/*.bigBed)
do
	# convert bigbed to bed
	bed=$(echo $bigbed | sed 's/.bigBed/.bed/')
	ATAC_seq/utils/bigBedToBed $bigbed $bed

	# merge overlapping regions in bed
	bed_merged=$(echo $bed | sed 's/.bed/.merged.bed/')
	cat $bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3-1,$4}' > ATAC_seq/data/bed.tmp
	bedtools merge -i ATAC_seq/data/bed.tmp -c 4 -o mean > ATAC_seq/data/bed.merged.tmp
	cat ATAC_seq/data/bed.merged.tmp | awk -F '\t' -v OFS='\t' '{print $1,$2,$3+1,$4}' > $bed_merged

	# remove temp files
	rm ATAC_seq/data/bed.tmp ATAC_seq/data/bed.merged.tmp

	# convert bed to bedgraph
	bedgraph=$(echo $bed | sed 's/.bed/.bedgraph/')
	cut -f1,2,3,4 $bed_merged | sort -k1,1 -k2,2n > $bedgraph

	# convert bedgraph to bigwig
	bigwig=$(echo $bedgraph | sed 's/.bedgraph/.bigWig/')
	ATAC_seq/utils/bedGraphToBigWig $bedgraph ATAC_seq/utils/mm10.chrom.sizes.txt $bigwig
done
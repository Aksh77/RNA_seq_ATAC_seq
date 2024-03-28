# Description: Preprocess the data for the analysis

for bigbed in $(ls atacseq_data/*/*.bigBed); do
	# convert bigbed to bed
	bed=$(echo $bigbed | sed 's/.bigBed/.bed/')
	utils/bigBedToBed $bigbed $bed

	# convert bed to bedgraph
	bedgraph=$(echo $bed | sed 's/.bed/.bedgraph/')
	cut -f1,2,3,4 $bed | sort -k1,1 -k2,2n > $bedgraph

	# merge overlapping regions in bedgraph
	bedgraph_merged=$(echo $bedgraph | sed 's/.bedgraph/.merged.bedgraph/')
	bedtools merge -i $bedgraph -c 4 -d 0 -o mean > $bedgraph_merged
	# remove the original bedgraph
	rm $bedgraph

	# convert bedgraph to bigwig
	bigwig=$(echo $bedgraph_merged | sed 's/.merged.bedgraph/.bigWig/')
	utils/bedGraphToBigWig $bedgraph_merged utils/mm10.chrom.sizes.txt $bigwig
done


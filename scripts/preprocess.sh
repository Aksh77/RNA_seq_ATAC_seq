# Description: Preprocess the data for the analysis

for bigbed in $(ls data/input_data/*/*.bigBed); do
	# convert bigbed to bed
	bed=$(echo $bigbed | sed 's/.bigBed/.bed/')
	utils/bigBedToBed $bigbed $bed

	# convert bed to bedgraph
	bedgraph=$(echo $bed | sed 's/.bed/.bedgraph/')
	cut -f1,2,3,4 $bed | sort -k1,1 -k2,2n > $bedgraph

	# merge overlapping regions in bedgraph
	bedgraph_merged=$(echo $bedgraph | sed 's/.bedgraph/.merged.bedgraph/')
	cat $bedgraph | awk -F '\t' -v OFS='\t' '{print $1,$2,$3-1,$4}' > bedGraph.tmp
	bedtools merge -i bedGraph.tmp -c 4 -o mean > bedGraph.merged.tmp
	cat bedGraph.merged.tmp | awk -F '\t' -v OFS='\t' '{print $1,$2,$3+1,$4}' > $bedgraph_merged
	
	# remove temp files
	rm bedGraph.tmp bedGraph.merged.tmp

	# convert bedgraph to bigwig
	bigwig=$(echo $bedgraph_merged | sed 's/.merged.bedgraph/.bigWig/')
	utils/bedGraphToBigWig $bedgraph_merged utils/mm10.chrom.sizes.txt $bigwig
done


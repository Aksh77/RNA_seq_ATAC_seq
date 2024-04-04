# # preprocess the data for the analysis
# for bigbed in $(ls data/input_data/*/*.bigBed); do
# 	# convert bigbed to bed
# 	bed=$(echo $bigbed | sed 's/.bigBed/.bed/')
# 	utils/bigBedToBed $bigbed $bed

# 	# merge overlapping regions in bed
# 	bed_merged=$(echo $bed | sed 's/.bed/.merged.bed/')
# 	cat $bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3-1,$4}' > data/bed.tmp
# 	bedtools merge -i data/bed.tmp -c 4 -o mean > data/bed.merged.tmp
# 	cat data/bed.merged.tmp | awk -F '\t' -v OFS='\t' '{print $1,$2,$3+1,$4}' > $bed_merged

# 	# remove temp files
# 	rm data/bed.tmp data/bed.merged.tmp
# done

# get consensus peaks from replicates for each cell line
for dir in $(ls -d data/input_data/*); do
	# get cell line name from subdirectory name
	seq=$(echo $dir | sed 's/data\/input_data\///')
	seq=$(echo $seq | sed 's/\///')
	seq=$(echo $seq | cut -d'_' -f2)
	echo $seq

	# output directory
	out_dir=data/consensus_peaks
	mkdir -p $out_dir

	# intersect bed files in the same subdirectory to get consensus peaks
	bed1=$(ls $dir/*.merged.bed | head -n 1)
	bed2=$(ls $dir/*.merged.bed | tail -n 1)
	consensus_bed=$out_dir/${seq}_consensus_peaks.bed
	bedtools intersect -a $bed1 -b $bed2  -f 0.50 -r > $consensus_bed

	# convert consensus bed to bedgraph
	bedgraph=$(echo $consensus_bed | sed 's/.bed/.bedgraph/')
	cut -f1,2,3,4 $consensus_bed | sort -k1,1 -k2,2n > $bedgraph

	# convert bedgraph to bigwig
	bigwig=$(echo $bedgraph | sed 's/.bedgraph/.bigWig/')
	utils/bedGraphToBigWig $bedgraph utils/mm10.chrom.sizes.txt $bigwig

	# get average signal for each peak for downstream analysis
	bed_peak_list=$(echo $bigwig | sed 's/.bigWig/.peak_list.bed/')
	tab_file=$(echo $bigwig | sed 's/.bigWig/.tab/')
	cat $consensus_bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > $bed_peak_list
	echo "utils/bigWigAverageOverBed $bigwig $bed_peak_list $tab_file"
	utils/bigWigAverageOverBed $bigwig $bed_peak_list $tab_file
done




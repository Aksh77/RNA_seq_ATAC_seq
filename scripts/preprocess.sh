# preprocess the data for the analysis
for bigbed in $(ls data/input_data/*/*.bigBed); do
	# convert bigbed to bed
	bed=$(echo $bigbed | sed 's/.bigBed/.bed/')
	utils/bigBedToBed $bigbed $bed

	# merge overlapping regions in bed
	bed_merged=$(echo $bed | sed 's/.bed/.merged.bed/')
	cat $bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3-1,$4}' > data/bed.tmp
	bedtools merge -i data/bed.tmp -c 4 -o mean > data/bed.merged.tmp
	cat data/bed.merged.tmp | awk -F '\t' -v OFS='\t' '{print $1,$2,$3+1,$4}' > $bed_merged

	# remove temp files
	rm data/bed.tmp data/bed.merged.tmp

	# convert bed to bedgraph
	bedgraph=$(echo $bed_merged | sed 's/.merged.bed/.bedgraph/')
	cut -f1,2,3,4 $bed_merged | sort -k1,1 -k2,2n > $bedgraph

	# convert bedgraph to bigwig
	bigwig=$(echo $bedgraph | sed 's/.bedgraph/.bigWig/')
	utils/bedGraphToBigWig $bedgraph utils/mm10.chrom.sizes.txt $bigwig

	# create output file path for tab file
	out_file=$(echo $bigwig | sed 's/data\/input_data/data\/output_data/')
	out_file=$(echo $out_file | sed 's/.bigWig/.tab/')
	echo $out_file
	mkdir -p $(dirname $out_file)

	# get average signal for each peak for downstream analysis
	bed_peak_list=$(echo $bed_merged | sed 's/.merged.bed/.peak_list.bed/')	
	tab_file=$(echo $bigwig | sed 's/.merged.bed/.tab/')
	cat $bed_merged | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > $bed_peak_list
	utils/bigWigAverageOverBed $bigwig $bed_peak_list $out_file
done

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
	bed1=$(ls $dir/*.bed | head -n 1)
	bed2=$(ls $dir/*.bed | tail -n 1)
	bedtools intersect -a $bed1 -b $bed2  -f 0.50 -r > $out_dir/${seq}_consensus_peaks.bed
done

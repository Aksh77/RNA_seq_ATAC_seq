# peakcalling using MACS2
for dir in $(ls -d ATAC_seq/data/input_data/*); do
	# get cell line name from subdirectory name
	seq=$(echo $dir | sed 's/ATAC_seq/data\/input_data\///')
	seq=$(echo $seq | sed 's/\///')
	seq=$(echo $seq | cut -d'_' -f2)
	echo $seq

	# call peaks using MACS2
	rep1=$(ls $dir/*.bam | head -n 1)
	rep2=$(ls $dir/*.bam | tail -n 1)
	peak_dir=ATAC_seq/data/peaks/$seq
	mkdir -p $peak_dir
	macs2 callpeak -t $rep1 -f BAM -g mm -n ${seq}_1 --outdir $peak_dir
	macs2 callpeak -t $rep2 -f BAM -g mm -n ${seq}_2 --outdir $peak_dir
done

# get consensus peaks from results of the replicates
SEQS="hsc cmp cfue erythroblast"
for seq in $SEQS; do
	seq=hsc
	dir=ATAC_seq/data/peaks/${seq}
	rep1=${dir}/${seq}_1_peaks.narrowPeak
	rep2=${dir}/${seq}_2_peaks.narrowPeak
	consensus=${dir}/${seq}_peaks.narrowPeak
	bedtools intersect -a $rep1 -b $rep2 -f 0.5 -r > $consensus
done



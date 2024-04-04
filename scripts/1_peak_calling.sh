# peakcalling using MACS2
for dir in $(ls -d data/input_data/*); do
	# get cell line name from subdirectory name
	seq=$(echo $dir | sed 's/data\/input_data\///')
	seq=$(echo $seq | sed 's/\///')
	seq=$(echo $seq | cut -d'_' -f2)
	echo $seq

	# call peaks using MACS2
	rep1=$(ls $dir/*.bam | head -n 1)
	rep2=$(ls $dir/*.bam | tail -n 1)
	peak_dir=data/peaks/$seq
	mkdir -p $peak_dir
	macs2 callpeak -t $rep1 $rep2 -f BAM -g mm -n $seq --outdir $peak_dir
done

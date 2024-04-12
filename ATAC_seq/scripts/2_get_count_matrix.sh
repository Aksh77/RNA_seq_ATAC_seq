# merge consensus peaks into peakset
hsc=data/peaks/hsc/hsc_peaks.narrowPeak
cfue=data/peaks/cfue/cfue_peaks.narrowPeak
hsc_cfue=data/peaks/hsc_cfue_peaks.narrowPeak
bedops -m $hsc $cfue > $hsc_cfue 

# get files for featureCounts
output_dir=data/output_data/featureCounts/hsc_cfue
mkdir -p $output_dir
hsc_cfue_bed=${output_dir}/hsc_cfue_peaks.bed
hsc_cfue_saf=${output_dir}/hsc_cfue_peaks.saf
hsc_cfue_counts=${output_dir}/hsc_cfue_peaks.counts
hsc_cfue_counts_tsv=${output_dir}/hsc_cfue_peaks.tsv
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="hsc_cfue_merged_"++nr;  print peakid,$1,$2,$3,"."}' $hsc_cfue > $hsc_cfue_saf
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="hsc_cfue_merged_"++nr;  print $1,$2,$3,peakid,"0","."}' $hsc_cfue > $hsc_cfue_bed

BAM1=data/input_data/1_hsc/ENCFF255IVU.bam
BAM2=data/input_data/1_hsc/ENCFF662DYG.bam
BAM3=data/input_data/3_cfue/ENCFF599ZDJ.bam
BAM4=data/input_data/3_cfue/ENCFF796ZSB.bam

# get counts matrix for peaks
featureCounts -F SAF -a $hsc_cfue_saf --fracOverlap 0.2 -o $hsc_cfue_counts $BAM1 $BAM2 $BAM3 $BAM4

# remove first line starting with #
awk '(NR>1)' $hsc_cfue_counts > $hsc_cfue_counts_tsv

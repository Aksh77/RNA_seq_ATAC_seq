# merge consensus peaks into peakset
hsc=data/consensus_peaks/hsc_consensus_peaks.bed
cfue=data/consensus_peaks/cfue_consensus_peaks.bed
hsc_cfue=data/consensus_peaks/hsc_cfue_merged_consensus_peaks.bed
bedops -m $hsc $cfue > $hsc_cfue 

# get files for downstream analysis
hsc_cfue_bed=data/consensus_peaks/hsc_cfue_peaks.bed
hsc_cfue_saf=data/consensus_peaks/hsc_cfue_peaks.saf
hsc_cfue_counts=data/output_data/hsc_cfue_peaks.counts
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="hsc_cfue_merged_"++nr;  print peakid,$1,$2,$3,"."}' $hsc_cfue > $hsc_cfue_saf
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="hsc_cfue_merged_"++nr;  print $1,$2,$3,peakid,"0","."}' $hsc_cfue > $hsc_cfue_bed

featureCounts -F SAF -a $hsc_cfue_saf --fracOverlap 0.2 -o $hsc_cfue_counts data/input_data/1_hsc/ENCFF255IVU.bam data/input_data/1_hsc/ENCFF662DYG.bam data/input_data/3_cfue/ENCFF599ZDJ.bam data/input_data/3_cfue/ENCFF796ZSB.bam

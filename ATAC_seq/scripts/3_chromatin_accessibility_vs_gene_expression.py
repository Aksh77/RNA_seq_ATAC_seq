import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# get contrast from command line argument
parser = argparse.ArgumentParser(description='Plot the relationship between chromatin accessibility and gene expression')
parser.add_argument('--contrast', type=str, help='Specify contrast in the format celltype1_celltype2, for example, hsc_cfue')
contrast = parser.parse_args().contrast

# get paths to differential accessibility, peak annotation, and gene expression data
DIFF_ACC = f"ATAC_seq/data/output_data/differential_accessibility/{contrast}/{contrast}_diff_acc.csv"
PEAK_ANNOT = f"ATAC_seq/data/output_data/peak_annotation/{contrast}/{contrast}_annotated_peaks.tsv"
GENE_EXP = f"ATAC_seq/data/output_data/gene_expression/{contrast}_gene_exp.csv"

# read differential chromatin accessibility data
diff_acc = pd.read_csv(DIFF_ACC, sep=",")
diff_acc["chr_start_end"] = diff_acc["seqnames"] + "_" + diff_acc["start"].astype(str) + "_" + diff_acc["end"].astype(str)
diff_acc.rename(columns={"FC": "Fold"}, inplace=True)
diff_acc.set_index("chr_start_end", inplace=True)
diff_acc = diff_acc[["Fold"]]

# get genes near differentially accessible peaks
peak_annot = pd.read_csv(PEAK_ANNOT, sep="\t")
peak_annot["chr_start_end"] = peak_annot["seqnames"] + "_" + peak_annot["start"].astype(str) + "_" + peak_annot["end"].astype(str)
peak_annot.set_index("chr_start_end", inplace=True)
peak_annot = peak_annot[["ENSEMBL", "SYMBOL", "annotation", "distanceToTSS", "annotation"]]

# join differential accessibility data with peak annotation data
diff_acc_gene = diff_acc.join(peak_annot)
diff_acc_gene.reset_index(inplace=True)
diff_acc_gene.sort_values("ENSEMBL", inplace=True)

# get gene expression data
gene_exp = pd.read_csv(GENE_EXP)
gene_exp.rename(columns={
    "Unnamed: 0": "gene_id",
    "log2FoldChange": "gene_exp_log2FC"
}, inplace=True)
gene_exp["gene_id"] = gene_exp["gene_id"].str.split(".").str[0]
gene_exp = gene_exp[["gene_id", "gene_exp_log2FC"]]

# join gene expression data with differential accessibility data on ENSEMBL and gene ID
diff_acc_gene_exp = diff_acc_gene.merge(gene_exp, left_on="ENSEMBL", right_on="gene_id", how='left')

# perform correlation analysis between differential accessibility and gene expression
corr = diff_acc_gene_exp[["Fold", "gene_exp_log2FC"]].corr()["Fold"]["gene_exp_log2FC"]

# plot the relationship between differential accessibility and gene expression
sns.set(style="whitegrid")
plt.figure(figsize=(7, 7))
# scatter plot
sns.scatterplot(
    x="Fold", 
    y="gene_exp_log2FC", 
    data=diff_acc_gene_exp,
    alpha=0.5
)
# fit a regression line
sns.regplot(
    x="Fold", 
    y="gene_exp_log2FC", 
    data=diff_acc_gene_exp,
    scatter=False,
    color='black'
)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("Differential Accessibility Fold change", fontsize=14, labelpad=10)
plt.ylabel("Gene Expression log2FC", fontsize=14, labelpad=10)
plt.title(f"Correlation coefficient: {corr:.2f}", fontsize=14, pad=10)
plt.tight_layout()
plt.savefig(f"ATAC_seq/data/output_data/chromatin_accessibility_vs_gene_expression/{contrast}.png")

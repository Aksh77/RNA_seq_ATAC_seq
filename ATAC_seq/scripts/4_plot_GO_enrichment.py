import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# get contrast from command line argument
parser = argparse.ArgumentParser(description='Gene Ontology enrichment plots for ATAC-seq data')
parser.add_argument('--contrast', type=str, help='Specify contrast in the format celltype1_celltype2, for example, hsc_cfue')
contrast = parser.parse_args().contrast

# plot top 20 molecular function GO terms
GO_MF = f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_MF.tsv"
go_mf = pd.read_csv(GO_MF, sep="\t")
go_mf = go_mf.sort_values('p.adjust', ascending=True)
go_mf = go_mf.head(20)

plt.figure(figsize=(4,6))
sns.barplot(data=go_mf, x='Count', y='Description', palette='viridis')
plt.ylabel('')
plt.yticks(fontsize=12)
plt.title(f'Top 20 enriched GO terms for {contrast} (Molecular Function)', pad=20)
plt.tight_layout()
plt.savefig(f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_MF.png", bbox_inches='tight')

# plot top 20 biological process GO terms
GO_BP = f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_BP.tsv"
go_bp = pd.read_csv(GO_BP, sep="\t")
go_bp = go_bp.sort_values('p.adjust', ascending=True)
go_bp = go_bp.head(20)

plt.figure(figsize=(4,6))
sns.barplot(data=go_bp, x='Count', y='Description', palette='viridis')
plt.ylabel('')
plt.yticks(fontsize=12)
plt.title(f'Top 20 enriched GO terms for {contrast} (Biological Process)', pad=20)
plt.tight_layout()
plt.savefig(f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_BP.png", bbox_inches='tight')
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# get contrast from command line argument
parser = argparse.ArgumentParser(description='Gene Ontology enrichment plots for ATAC-seq data')
parser.add_argument('--contrast', type=str, help='Specify contrast in the format celltype1_celltype2, for example, hsc_cfue')
contrast = parser.parse_args().contrast

# plot top 10 molecular function GO terms 
GO_MF_UP = f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_MF_upregulated.tsv"
go_mf_upregulated = pd.read_csv(GO_MF_UP, sep="\t")
go_mf_upregulated = go_mf_upregulated.sort_values('p.adjust', ascending=True)
go_mf_upregulated = go_mf_upregulated.head(20)
plt.figure(figsize=(4,6))
sns.barplot(data=go_mf_upregulated, x='Count', y='Description', palette='viridis')
plt.ylabel('')
plt.yticks(fontsize=12)
plt.title(f'Top 20 enriched GO terms for {contrast} (Molecular Function)', pad=20)
plt.tight_layout()
plt.savefig(f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_MF_up.png", bbox_inches='tight')

GO_MF_DOWN = f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_MF_downregulated.tsv"
go_mf_downregulated = pd.read_csv(GO_MF_DOWN, sep="\t")
go_mf_downregulated = go_mf_downregulated.sort_values('p.adjust', ascending=True)
go_mf_downregulated = go_mf_downregulated.head(20)
plt.figure(figsize=(4,6))
sns.barplot(data=go_mf_downregulated, x='Count', y='Description', palette='viridis')
plt.ylabel('')
plt.yticks(fontsize=12)
plt.title(f'Top 20 enriched GO terms for {contrast} (Molecular Function)', pad=20)
plt.tight_layout()
plt.savefig(f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_MF_down.png", bbox_inches='tight')

# plot top 20 biological process GO terms
GO_BP_UP = f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_BP_upregulated.tsv"
go_bp_upregulated = pd.read_csv(GO_BP_UP, sep="\t")
go_bp_upregulated = go_bp_upregulated.sort_values('p.adjust', ascending=True)
go_bp_upregulated = go_bp_upregulated.head(20)
plt.figure(figsize=(4,6))
sns.barplot(data=go_bp_upregulated, x='Count', y='Description', palette='viridis')
plt.ylabel('')
plt.yticks(fontsize=12)
plt.title(f'Top 20 enriched GO terms for {contrast} (Biological Process)', pad=20)
plt.tight_layout()
plt.savefig(f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_BP_up.png", bbox_inches='tight')

GO_BP_DOWN = f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_BP_downregulated.tsv"
go_bp_downregulated = pd.read_csv(GO_BP_DOWN, sep="\t")
go_bp_downregulated = go_bp_downregulated.sort_values('p.adjust', ascending=True)
go_bp_downregulated = go_bp_downregulated.head(20)
plt.figure(figsize=(4,6))
sns.barplot(data=go_bp_downregulated, x='Count', y='Description', palette='viridis')
plt.ylabel('')
plt.yticks(fontsize=12)
plt.title(f'Top 20 enriched GO terms for {contrast} (Biological Process)', pad=20)
plt.tight_layout()
plt.savefig(f"ATAC_seq/data/output_data/GO_enrichment/{contrast}/{contrast}_enriched_GO_BP_down.png", bbox_inches='tight')
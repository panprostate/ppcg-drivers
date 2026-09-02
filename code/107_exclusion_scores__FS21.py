import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


genes = pd.read_csv('~/exclusion_score_RGAs.tsv', 
                     sep='\t', 
                    comment=None,
                    index_col=False,
                    header=0        
                    )

genes = genes.sort_values('ratio', ascending=False)
genes = genes[~genes.gene.str.contains('SV|SNV')]

df_sorted = genes.sort_values('ratio', ascending=False)

plt.figure(figsize=(6, 2))

bars = plt.bar(df_sorted["gene"], df_sorted["ratio"], color='#444444', linewidth=0.5, width=0.6) 

plt.axhline(1, color='black', linestyle='--', linewidth=0.7)

plt.ylabel("Exclusivity score ratio\n(Observed / Expected)", fontsize=6, fontname='Arial')
plt.xticks(rotation=90, fontsize=6, fontname='Arial', ha='right') 

ax = plt.gca()
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))
ax.tick_params(axis='y', labelsize=6)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

colors = ['black'] * len(df_sorted)

red_genes = df_sorted["gene"].isin(genes[genes.pvalue_higher < 0.05]["gene"])
blue_genes = df_sorted["gene"].isin(genes[genes.pvalue_lower < 0.05]["gene"])

for i, (bar, gene) in enumerate(zip(bars, df_sorted["gene"])):
    if red_genes.iloc[i]:
        colors[i] = 'r'
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05, '*', 
                ha='center', va='top', fontsize=8, fontweight='bold', color='r')
    elif blue_genes.iloc[i]:
        colors[i] = 'blue'
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05, '*', 
                ha='center', va='top', fontsize=8, fontweight='bold', color='blue')

xticks = ax.get_xticklabels()
for i, label in enumerate(xticks):
    label.set_color(colors[i])

plt.xlim(-1, len(df_sorted["gene"]) + 0)

#plt.savefig("~/exclusion_ratio.png", format="png", bbox_inches='tight', dpi=1200)

plt.show()
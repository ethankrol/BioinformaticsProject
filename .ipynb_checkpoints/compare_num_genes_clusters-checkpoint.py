import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt

num_genes_list = [10, 100, 1000, 5000, 10000]

###################################################
## Chi-Squared Results
###################################################

comparison_results = [[None] * (len(num_genes_list)+1) for _ in range(len(num_genes_list)+1)]

for i, num_genes in enumerate(num_genes_list):
    comparison_results[i+1][0] = num_genes
    comparison_results[0][i+1] = num_genes

comparison_results[0][0] = 'Num Genes'


# sample each pair from num_genes_list
for i in range(len(num_genes_list)):
    for j in range(i+1, len(num_genes_list)):
        num_genes1 = num_genes_list[i]
        num_genes2 = num_genes_list[j]
        
        # read the data
        df1 = pd.read_csv(f"results_assignment3/2_spectral_clusters_from_{num_genes1}_genes.tsv", sep='\t')
        df2 = pd.read_csv(f"results_assignment3/2_spectral_clusters_from_{num_genes2}_genes.tsv", sep='\t')

        # contingency table
        contingency_table = pd.crosstab(df1['Cluster'], df2['Cluster'])

        # Perform the Chi-Square test of independence
        chi2, p_value, dof, expected = chi2_contingency(contingency_table)

        p_value = f'{float(p_value):.3e}'

        # add results to table
        comparison_results[i+1][j+1] = p_value


# Create a figure and axis
fig, ax = plt.subplots()

# Hide axes
ax.axis('tight')
ax.axis('off')

# Create a table and add it to the axis
table = ax.table(cellText=comparison_results, loc='center', cellLoc='center')

table.auto_set_font_size(False)
table.set_fontsize(13)
table.scale(1.5, 1.5) 

plt.title('Comparison of Spectral Clustering Results\nwith Different Numbers of Genes Used in Clustering', fontsize=16, y=0.8)

caption = "This table displays the p-values obtained from pairwise Chi-Squared tests."
plt.text(0, -0.03, caption, ha='center', va='top', fontsize=12, wrap=True)

# Save the table as a PNG file
plt.savefig('results_assignment3/spectral_num_genes_comparison.png', bbox_inches='tight', dpi=300)

###################################################
## Sankey Diagram
###################################################

import plotly.graph_objects as go

source = []
target = []
value = []

for i in range(len(num_genes_list) - 1):
    num_genes1 = num_genes_list[i]
    num_genes2 = num_genes_list[i+1]
    
    # read the data
    df1 = pd.read_csv(f"results_assignment3/2_spectral_clusters_from_{num_genes1}_genes.tsv", sep='\t')
    df2 = pd.read_csv(f"results_assignment3/2_spectral_clusters_from_{num_genes2}_genes.tsv", sep='\t')

    # contingency table
    contingency_table = pd.crosstab(df1['Cluster'], df2['Cluster'])

    for j in range(contingency_table.shape[0]):
        for k in range(contingency_table.shape[1]):
            source.append(2 * i + j)
            target.append(2 * i + 2 + k)
            value.append(int(contingency_table.iloc[j, k]))

# Create a Sankey diagram
fig = go.Figure(data=[go.Sankey(
node = dict(
    pad = 50,
    thickness = 10,
    line = dict(color = "black", width = 0.5),
    label = [f"{num_genes} genes<br>Cluster {i}" for num_genes in num_genes_list for i in range (1,3)],
    color = "blue"
),
link = dict(
    source = source,
    target = target,
    value = value,
))])

fig.update_layout(title_text="Sankey Diagram of Spectral Clustering Results\nwith Different Numbers of Genes Used in Clustering", font_size=10)
# save the figure as png
fig.write_image("results_assignment3/spectral_num_genes_sankey.png")

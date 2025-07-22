"""
Simulates how the number of contributing genes increases the variance in disease liability (PGS). 
Generates a 2-panel figure illustrating these relationships.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

np.random.seed(0)

# Function to calculate PGS (additive + epistatic) for an individual
def calculate_PGS(SNP_AF, SNP_eff, ploidy, interactions):
    genotypes = np.random.binomial(ploidy, SNP_AF)
    additive_score = np.dot(SNP_eff, genotypes)
    epistatic_score = np.einsum('i,ij,j->', genotypes, np.triu(interactions), genotypes)
    return additive_score + epistatic_score, genotypes, additive_score

# Parameters
nSNPs = [x * 50 + 50 for x in range(40)]
coloured_snps = [500, 1000, 2000]
nPop = 1000
nReplicates = 10

variance = []
variance_epi = []
nSNPs_for_plot = []
PGS_replicates = []

for n in range(len(nSNPs)):
    print(n)
    for rep in range(nReplicates):
        SNP_AF = np.random.uniform(0.0, 1.0, nSNPs[n])
        SNP_eff = np.random.normal(0.0, 0.01, nSNPs[n])
        epistatic_interactions = np.random.normal(0.0, 0.0004, (nSNPs[n], nSNPs[n]))

        PGS_population = []
        PGSepi_population = []
        for i in range(nPop):
            PGS_epi, genotype, PGS = calculate_PGS(SNP_AF, SNP_eff, ploidy=2, interactions=epistatic_interactions)
            PGSepi_population.append(PGS_epi)
            PGS_population.append(PGS)

        if nSNPs[n] in coloured_snps:
            PGS_replicates.append(PGS_population)

        variance.append(np.var(PGS_population))
        variance_epi.append(np.var(PGSepi_population))
        nSNPs_for_plot.append(nSNPs[n])

# Plotting style specifications
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=SMALL_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# --- Left Plot (KDE for disease liability) ---
ax = axes[0]
for rep in range(nReplicates):
    rescaled = PGS_replicates[rep] - np.mean(PGS_replicates[rep])
    rescaled2 = PGS_replicates[rep + nReplicates] - np.mean(PGS_replicates[rep + nReplicates])
    rescaled3 = PGS_replicates[rep + 2 * nReplicates] - np.mean(PGS_replicates[rep + 2 * nReplicates])
    sns.kdeplot(rescaled, alpha=0.5, color='green', ax=ax)
    sns.kdeplot(rescaled2, alpha=0.5, color='blue', ax=ax)
    sns.kdeplot(rescaled3, alpha=0.5, color='red', ax=ax)

# Gradient background
colors = ["white", "#DD908F"]
cmap = LinearSegmentedColormap.from_list("custom_gradient", colors)
gradient = np.hstack([
    np.linspace(0, 0, 100),
    np.linspace(0, 1, 100)
]).reshape(1, -1)
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
ax.imshow(gradient, extent=[xmin, xmax, ymin, ymax], aspect='auto', cmap=cmap, alpha=1.0, zorder=1)

# Legend
line1 = plt.Line2D([], [], color='green')
line2 = plt.Line2D([], [], color='blue')
line3 = plt.Line2D([], [], color='red')
custom_lines = [line1, line2, line3]
ax.legend(custom_lines, ['500', '1000', '2000'], loc='upper right', title='Number of genes')

ax.set_xlabel('Additive disease liability \n (Polygenic score)')
ax.set_ylabel('')
ax.set_yticks([])

# --- Right Plot: increasing variance of disease liability ---
ax = axes[1]

# Color and style setup
color_by_type = {"Additive": "black", "Additive and epistatic": "#CC6600"}
highlight_colors = {500: 'green', 1000: 'blue', 2000: 'red'}

# Prepare long-form variance data
records = []
for i in range(len(nSNPs_for_plot)):
    records.append({
        "nSNPs": nSNPs_for_plot[i],
        "Variance": variance[i],
        "Type": "Additive"
    })
    records.append({
        "nSNPs": nSNPs_for_plot[i],
        "Variance": variance_epi[i],
        "Type": "Additive and epistatic"
    })
df = pd.DataFrame(records)

# Group by type and SNP count
summary = (
    df.groupby(["nSNPs", "Type"])
      .agg(mean_var=("Variance", "mean"),
           se_var=("Variance", lambda x: np.std(x, ddof=1) / np.sqrt(len(x))))
      .reset_index()
)

# Plot mean + shaded SE for each effect type
for effect_type in ["Additive", "Additive and epistatic"]:
    subset = summary[summary["Type"] == effect_type]
    ax.plot(
        subset["nSNPs"],
        subset["mean_var"],
        color=color_by_type[effect_type],
        label=effect_type,
        linewidth=2
    )
    ax.fill_between(
        subset["nSNPs"],
        subset["mean_var"] - subset["se_var"],
        subset["mean_var"] + subset["se_var"],
        color=color_by_type[effect_type],
        alpha=0.3
    )

for i in range(len(nSNPs_for_plot)):
    if nSNPs_for_plot[i] in coloured_snps:
        ax.scatter(nSNPs_for_plot[i], variance[i], color=highlight_colors[nSNPs_for_plot[i]], alpha = 0.5, zorder = 3)

# Labels and legend
ax.set_xlabel("Number of genes")
ax.set_ylabel("Disease liability variance")
ax.legend(title="Interaction type", loc="upper left")

ax.set_xlim(left=0.0)
ax.set_ylim(bottom=0.0)

# Labels
plt.gcf().text(0.0, 1.0, 'A)', fontsize=14, fontweight='bold', va='top', ha='left')
plt.gcf().text(0.5, 1.0, 'B)', fontsize=14, fontweight='bold', va='top', ha='left')

plt.tight_layout()
plt.show()
plt.savefig("pgs-distributions.png", dpi=300, bbox_inches="tight")
plt.close()


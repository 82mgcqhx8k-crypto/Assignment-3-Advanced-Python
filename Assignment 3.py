# Assignment 3: Advanced Python
# Yash
# 10/05/26

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1a. Load expression, gene annotation, and sample metadata
df_expression = pd.read_excel("Gene_Expression_Data.xlsx", index_col = 0)
df_gene = pd.read_csv("Gene_Information.csv")
df_sample = pd.read_csv("Sample_Information.tsv", sep = "\t")

# Map each sample ID to its phenotype group
change_sample_name = df_sample["group"].to_dict()

# Rename expression columns with phenotype labels and sample numbers
new_columns = []
counts = {}

for col in df_expression.columns:
    pheno = str(change_sample_name.get(col, "unknown")).lower()
    counts[pheno] = counts.get(pheno, 0) + 1
    new_columns.append(f"{pheno}_{counts[pheno]}")

df_expression.columns = new_columns

# Separate tumor and normal samples
tumor_column = [c for c in df_expression.columns if c.startswith("tumor")]
normal_column = [c for c in df_expression.columns if c.startswith("normal")]

df_tumor = df_expression[tumor_column]
df_normal = df_expression[normal_column]


# Calculate each probe's average expression in each group
mean_tumor = df_tumor.mean(axis = 1)
mean_normal = df_normal.mean(axis = 1)

# Calculate relative expression change between tumor and normal samples
fold_change = (mean_tumor - mean_normal) / mean_normal

# Combine group means and fold change into one analysis table
df_analysis = pd.DataFrame({
    "Mean_Tumor": mean_tumor,
    "Mean_Normal": mean_normal,
    "Fold_Change": fold_change
})

# Make the probe IDs a regular column for merging with gene annotations
df_analysis = df_analysis.reset_index()
probe_column = df_analysis.columns[0]

# Add gene annotation information to the expression results
df_merged = pd.merge(
    df_analysis,
    df_gene,
    left_on = probe_column,
    right_on = "Probe_ID"
)

# Keep probes with an absolute fold change greater than 5
df_filtered = df_merged[df_merged["Fold_Change"].abs() > 5].copy()

# Label the group with the higher average expression for each probe
df_filtered["Higher_Expressed_In"] = np.where(
    df_filtered["Mean_Tumor"] > df_filtered["Mean_Normal"],
    "Tumor",
    "Normal"
)

# Retrieve expression values for the filtered differentially expressed probes
deg_probe = df_filtered["Probe_ID"].values
df_raw_filtered = df_expression.loc[deg_probe]

# Print results, line breaker to separate the two answers
print("  DEG Statistical Summary ")
if not df_filtered.empty:
    print(df_filtered.describe(include = "all"),'\n')
    print("\n--- Structural Regulation Balance ---")
    print(df_filtered["Higher_Expressed_In"].value_counts())

# Sort chromosome labels numerically where possible, then alphabetically
df_filtered["Chromosome"] = (
    df_filtered["Chromosome"].astype(str).str.strip()
)

def sort_key(chrom):
    # Remove a 'chr' prefix before sorting
    c = chrom.lower().replace("chr", "")
    if c.isdigit():
        return 0, int(c)
    return 1, c

sorted_chrom = sorted(
    df_filtered["Chromosome"].unique(),
    key = sort_key
)

# Plot the number of filtered DEGs on each chromosome
plt.figure(figsize=(14, 6))
sns.countplot(
    data = df_filtered,
    x = "Chromosome",
    order = sorted_chrom,
    palette = "Reds_r"
)
plt.title(
    "Distribution of DEGs by Chromosome",
    fontsize = 14,
    fontweight = "bold"
)
plt.xlabel("Chromosome", fontsize=12)
plt.ylabel("# OF DEGs", fontsize=12)
plt.xticks(rotation=45)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()
plt.close()

# Plot chromosome counts split by the group with higher expression
fig, ax = plt.subplots(figsize=(14, 6.5))

# Assign a distinct color to each higher-expression group
group_colors = {"Tumor": "Red", "Normal": "Blue"}

sns.countplot(
    data = df_filtered,
    x = "Chromosome",
    hue = "Higher_Expressed_In",
    order = sorted_chrom,
    palette = group_colors,
    edgecolor = "white",
    linewidth = 0.8,
    alpha = 0.9,
    ax = ax
)
ax.set_title(
    "Genomic Distribution of DEGs by Chromosome",
    fontsize = 15,
    fontweight = "bold",
    pad = 15
)
ax.set_xlabel("Chromosome Location", fontsize = 12, labelpad = 10)
ax.set_ylabel("# ofFiltered DEGs", fontsize = 12, labelpad = 10)
plt.tight_layout()
plt.show()
plt.close()


# Calculate the percentage of DEGs with higher expression in each group
plt.figure(figsize = (6, 5))
reg_count = (
    df_filtered["Higher_Expressed_In"].value_counts(normalize = True) * 100
)

# tumor-higher probes = upregulated/normal-higher probes = downregulated
labels = [
    "Upregulated\n(Tumor > Normal)",
    "Downregulated\n(Normal > Tumor)"
]
percentages = [reg_count.get("Tumor", 0), reg_count.get("Normal", 0)]

# Plot group percentages and label each bar with its value
bars = plt.bar(
    labels,
    percentages,
    color = ["Red", "Blue"],
    edgecolor = "black",
    width = 0.5
)
plt.title(
    "Percentage of DEGs Regulated in Tumor Samples",
    fontsize = 13,
    fontweight = "bold"
)
plt.ylabel("Percentage of Total DEGs (%)", fontsize=11)
plt.ylim(0, 110)

for bar in bars:
    yval = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        yval + 2,
        f"{yval:.1f}%",
        ha = "center",
        va = "bottom",
        fontweight = "bold"
    )

plt.tight_layout()
plt.show()
plt.close()


# Standardize each probe's expression values to z-scores for visualization
plt.figure(figsize = (10, 8))
df_zscore = df_raw_filtered.apply(
    lambda x: (x - x.mean()) / x.std(),
    axis = 1
)

# Display the standardized expression patterns as a heatmap
sns.heatmap(
    df_zscore,
    cmap = "RdBu_r",
    center = 0,
    cbar_kws = {"label": "Z-score normalized expression"},
    xticklabels = True,
    yticklabels = len(df_zscore) <= 50
)
plt.title(
    "Heatmap of DEG Expression Profiles",
    fontsize = 13,
    fontweight = "bold"
)
plt.xlabel("Samples", fontsize = 12)
plt.ylabel("Probes", fontsize = 12)
plt.tight_layout()
plt.show()
plt.close()

# Color sample labels by phenotype in the clustermap
col_colors = [
    "#d95f02" if "tumor" in c else "#7570b3"
    for c in df_zscore.columns
]

# Cluster samples and probes by their standardized expression profiles
cluster_map = sns.clustermap(
    df_zscore,
    cmap = "RdBu_r",
    center = 0,
    col_colors = col_colors,
    cbar_kws = {"label": "Z-score"},
    figsize = (11, 9)
)
cluster_map.figure.suptitle(
    "Clustermap of Samples and Top DEGs",
    y = 1.02,
    fontsize = 14,
    fontweight = "bold"
)
plt.show()
plt.close()
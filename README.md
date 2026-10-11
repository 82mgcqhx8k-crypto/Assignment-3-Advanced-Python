# Assignment-3-Advanced-Python

# Introduction
Programmer: Yash
Language: Python, Python libraries used within the python language

# Purpose
The purpose of this assignment was to analyze tumor and normal gene expression data(sample.tsv), calculates relative expression changes, identifies probes with an absolute fold change greater than 5, and produces chromosome distribution plots and expression heatmaps

# Files required
First make sure to place these files in the same folder as the directory of Assignment 3

Gene_Expression_Data.xlsx — contains expression measurements, with unique various Probe_IDs listed as well

Gene_Information.csv — includes various gene expression levels, what chromosome is active, and different cytoband along with Gene_ID  

Sample_Information.tsv — has sample metadata, it shows what patient has tumor and what patient is normal, also their corresponding GSM 

# Explanation of what each important section of code does
1. Loads and prepares the data, so the code can read the expression table, sample metadata, and probe annotation table. It matches sample IDs from the Sample_Information dataset

2. Compares tumor and normal expression, calculates the average expression for each probe in each group, then calculates the tumor to normal fold change as stated in the assignment tab. Probes with an absolute fold change greater than 5 are flagged as differentially expressed genes or also known as DEGs under this threshold

3. Adds gene annotations, joins the probe level results to Gene_Information using Probe_ID, so results can be interpreted using gene and chromosome information

4. Summarizes and visualizes the results, prints summary statistics and the number of probes passing the fold change threshold. It also creates various plots showing the chromosome distribution of flagged probes and heatmaps of expression patterns

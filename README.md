# CDR-Finder
This repository contains a snakemake to identify and quantify Centromere Demethylation Regions (CDR).

This is done by:
- Bin into 5kbp windows
- Intersect with ALR Annotations from Repeat Masker
- Only take bins that drop below the threshold (median frequency in the bins with ALR)
- Merge consecutive bins
- Filter for bins greater than 50kbp

# Input
- target_bed: Bed File of Target Region Coordinates
- meth_tsv: Nanopolish methylation frequency TSV
- rm_out: Repeat Masker .out file

# Output

Based on median methylation profile across region, idenitfies and bins regions with methylation frequency below median.
```
results/{sample}_CDR.bed
```

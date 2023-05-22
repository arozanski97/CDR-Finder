# CDR-Finder
This repository contains a snakemake to identify and quantify hypomethylated regions within centromeres, or Centromere Dip Regions (CDRs; Gershman et al., Science, 2022).

This is done by:
- Bin into 5 kbp windows
- Intersect with ALR annotations from RepeatMasker
- Only take bins that fall below the threshold (median frequency in the bins annotated as “ALR”)
- Merge consecutive bins
- Filter for bins greater than 50 kbp

# Input
- target_bed: BED file of target region coordinates
- meth_tsv: Nanopolish methylation frequency TSV
- rm_out: RepeatMasker .out file

# Output
Based on the median methylation frequency across the region, identifies and bins regions with methylation frequency below median spanning >50 kbp.
```
results/{sample}_CDR.bed
```

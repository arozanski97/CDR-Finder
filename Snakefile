import pandas as pd
import numpy as np
from Bio import SeqIO

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile: 'config.yaml'

TARGET_BED = config['target_bed']
METH_TSV = config['meth_tsv']
RM_out = config['rm_out']
WINDOW = config['window_size']
SAMPLE = config['sample']



rule all:
	input:
		expand("results/{sample}_CDR.bed", sample=SAMPLE)

rule calc_windows:
	input:
		methylation_tsv = METH_TSV,
		target_bed = TARGET_BED
	output:
		binned_freq = "results/{sample}_binned_freq.bed"
	resources: 
		mem = 8,
		hrs = 1
	threads: 1
	params:
		window_size = WINDOW
	run:
		meth_df = pd.read_csv(input.methylation_tsv, header=0, sep="\t")
		meth_df_sorted = meth_df.sort_values(by=['start','end'])
		keep_window = pd.DataFrame()
		window_size = int(params.window_size)
		with open(input.target_bed, "r") as infile:
			for line in infile:
				chrom, start, stop = line.split("\t")[0], line.split("\t")[1], line.split("\t")[2]
				meth_subset = (meth_df_sorted.loc[(meth_df_sorted['chromosome'] == chrom ) & (meth_df_sorted['start'] >= int(start)) & (meth_df_sorted['end'] <= int(stop))]).reset_index()
				window_avg = pd.DataFrame(columns=['Chr','Bin','Freq'])
				start_track = int(start)
				stop_track = int(start) + window_size
				track_avg = []	
				for idx in meth_subset.index:
					if int(meth_subset.at[idx,'start'])  >= stop_track:
						RM_start = int(start_track)-int(start)
						RM_stop = int(stop_track)-int(start)
						window_avg =window_avg.append({'Chr':str(chrom), 'Bin':str(RM_start)+"-"+str(RM_stop), 'Freq':np.mean(track_avg)}, ignore_index=True)
						start_track = stop_track
						stop_track = start_track+window_size
						track_avg = []
					avg = [float(meth_subset.at[idx, 'methylated_frequency'])] * int(meth_subset.at[idx,'num_motifs_in_group'])
					track_avg.extend(avg)

				window_avg = window_avg.round({'Freq':2})
				keep_window = keep_window.append(window_avg)

		keep_window.reset_index(inplace=True)
		with open(output.binned_freq, "w+") as outfile:
			for idx_2 in keep_window.index:
				
				pos1 = str(keep_window.at[idx_2,'Bin']).split("-")[0]
				pos2 = str(keep_window.at[idx_2,'Bin']).split("-")[1]
				freq = str(keep_window.at[idx_2,'Freq'])
				chrom = str(keep_window.at[idx_2,'Chr'])
				outfile.write(chrom+"\t"+pos1+"\t"+pos2+"\t"+str(freq)+"\t"+str(idx_2)+"\n")

rule format_RM:
	input:
		repeat_masker = RM_out,
	output:
		repeat_bed = "results/{sample}_rm.bed"
	resources:
		mem = 8,
		hrs = 1,
	threads: 1
	shell:
		'''
		awk -v OFS="\t" '{{print $5, $6, $7, $10, $9}}' {input.repeat_masker} > {output.repeat_bed}
		'''
rule filter_RM:
	input:
		repeat_masker = rules.format_RM.output.repeat_bed
	output:
		rm_bed = "results/{sample}_rm_ALR.bed"
	resources:
		mem = 8,
		hrs = 1
	threads: 1
	run:
		RM_df = pd.read_csv(input.repeat_masker, header=None, sep="\t", names=['chr','start','stop','ALR','or'])
		RM_df['chr']=RM_df['chr'].str.split(pat=":").str[0]
		RM_df.to_csv(output.rm_bed, header=None, index=None, sep="\t")

rule intersect_RM:
	input:
		repeat_masker = rules.filter_RM.output.rm_bed,
		binned_freq = rules.calc_windows.output.binned_freq
	output:
		intersect_bed = "results/{sample}_intersect.bed",
		merged = "results/{sample}_rm_merged.bed"
	resources:
		mem = 8,
		hrs = 1
	threads: 1
	shell:
		'''
		module load bedtools/2.29.2
		bedtools merge -i {input.repeat_masker} -d 500 > {output.merged}
		bedtools intersect -a {input.binned_freq} -b {output.merged} -f 1 -wa  -u > {output.intersect_bed}
		'''
rule in_threshold:
	input:
		intersect_bed = rules.intersect_RM.output.intersect_bed,
		target_bed = TARGET_BED
	output:
		final_call = "results/{sample}_CDR.bed"
	resources:
		mem = 8,
		hrs = 1
	threads: 1
	run:
		intersect_bed = pd.read_csv(input.intersect_bed, header=None, sep="\t", names=['chrom','start','stop','freq','idx'])
		
		with open(input.target_bed) as infile, open(output.final_call, "w+") as outfile:
			for line in infile:
				window_avg = pd.DataFrame(columns=['Bin','Freq'])
				chrom, start, stop = line.split("\t")[0], line.split("\t")[1], line.split("\t")[2]
				intersect_bed_sub = intersect_bed.loc[intersect_bed['chrom'] == chrom]
				for index in intersect_bed_sub.index:
					chrom_int = intersect_bed_sub.at[index,'chrom']
					start_int = int(intersect_bed_sub.at[index,'start']) +int(start)
					stop_int = int(intersect_bed_sub.at[index,'stop']) +int(start)
					freq = intersect_bed_sub.at[index,'freq']
					idx = intersect_bed_sub.at[index,'idx']
					window_avg =window_avg.append({'Bin':str(start_int)+"-"+str(stop_int), 'Freq':float(freq), 'idx_2': int(idx)}, ignore_index=True)
				threshold = np.median(window_avg['Freq'].tolist())
				print(threshold)
				window_avg = window_avg.set_index('idx_2')
				window_avg = window_avg.loc[window_avg['Freq'] <= threshold]
				s = pd.Series(window_avg.index.tolist())
				groups = s.groupby(s.diff().ne(1).cumsum()).apply(lambda x: [x.iloc[0], x.iloc[-1]] if len(x) >= 2 else [x.iloc[0]]).tolist()
				groups = [sub for sub in groups if len(sub)>1]
				for y in groups:
					pos1 = str(window_avg.at[y[0],'Bin']).split("-")[0]
					pos2 = str(window_avg.at[y[1],'Bin']).split("-")[1]
					if int(pos2)-int(pos1) > 50000:
						outfile.write(chrom+"\t"+pos1+"\t"+pos2+"\n")

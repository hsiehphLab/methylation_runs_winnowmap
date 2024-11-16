import pandas as pd
import os
import numpy as np
import pysam
import gzip
from Bio import SeqIO
import re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
#configfile: 'config.yaml'

#REF = config['REF']
MANIFEST = 'methyl.tsv'

manifest_df = pd.read_csv(MANIFEST, sep='\t', header=0)
manifest_df.set_index(['sample'], inplace=True, drop=False)

wildcard_constraints:


wildcard_constraints:
	hap = '|'.join(['hap1', 'hap2','unphased']),
	tech = '|'.join(['ont', 'hifi'])

tech_dict = {'ont': 'map-ont', 'hifi':'map-pb'}

rule all:
	input:
		[expand('results/{sample}/cpg/{sample}_{tech}_{hap}.bed', sample=SM, tech=manifest_df.at[SM, 'tech_type'], hap='unphased') for SM in (manifest_df.loc[manifest_df['unphased'] ==True].index)], 
		#[expand('results/{sample}/cpg/{tech}_{hap}.bed', sample=SM, hap=['hap1','hap2'], ) for SM in (manifest_df.loc[manifest_df['unphased'] == False].index)]

#GET LIST OF READS FOR EACH HAPLOTYPE
rule gethap:
	input:
		fasta = lambda wildcards: manifest_df.at[wildcards.sample,f'{wildcards.hap}_fasta']
	output:
		name_list = "results/{sample}/haplotype/{tech}_{hap}.list"
	resources:
		mem_mb = 12000,
		hrs = 12,
		threads = 1
	shell:
		'''
		zcat {input.fasta} |  grep -e ">" | sed -e 's/>//g' | cut -f1 -d " " > {output.name_list}
		'''
#FOR EACH HAP LOCATE READS IN BAM FILES
checkpoint find_bams:
	input:
		name_list = "results/{sample}/haplotype/{tech}_{hap}.list",
		fastq_fofn = lambda wildcards: manifest_df.at[wildcards.sample,'fastq_fofn'],
		bam_fofn = lambda wildcards: manifest_df.at[wildcards.sample,'bam_fofn'],
	output:
		readmanifest = touch("results/{sample}/haplotype/.{tech}_{hap}_split.done")
	resources:
		mem_mb = 12000,
		hrs = 12,
		threads = 1
	run:
		hap_reads = pd.read_csv(input.name_list,header=None, sep="\t", names=['read'])
		fastq_fofn = pd.read_csv(input.fastq_fofn, header=None, sep="\t", names=['file'])
		bam_fofn = pd.read_csv(input.bam_fofn, header=None, sep="\t", names=['file'])
		bam_fofn['movie'] = bam_fofn['file'].str.split("/").str[-1].str.replace(".bam",".fastq.gz")
		bam_fofn.set_index(['movie'], inplace=True, drop=False)
		file_df = pd.DataFrame()
		for fasta in fastq_fofn['file']:
			fai = fasta+".fai"
			fai_df = pd.read_csv(fai, header=None, names=['read','len','f1','f2','f3','f4'], sep="\t")
			fai_df['bam_path'] = bam_fofn.at[fasta.split("/")[-1], 'file']
			fai_df_hap = hap_reads.set_index('read').join(fai_df.set_index('read'))
			file_df = fai_df_hap.dropna()
			file_name = os.path.dirname(output.readmanifest)+"/"+wildcards.hap+"_"+fasta.split("/")[-1].replace('.fastq.gz','.list')
			file_df[['bam_path']].to_csv(file_name, header=False, sep="\t")

def find_reads(wildcards):
	if wildcards.hap == 'unphased':
		fofn = manifest_df.at[wildcards.sample,'bam_fofn']
		print(fofn)
		fofn_df = pd.read_csv(fofn, header=None, sep="\t", names=['file'])
		fill = range(len(fofn_df))
		flist = expand("results/{sample}/alignments/{tech}_{hap}_{name}_to_ref.bam", sample=wildcards.sample, hap=wildcards.hap, name=fill, tech=wildcards.tech)
		return {'bams':flist}
		#print(expand("results/{sample}/alignments/{tech}_{hap}_{name}_to_ref.bam", sample=wildcards.sample, hap=wildcards.hap, name=fill, tech=wildcards.tech))
		#return expand("results/{sample}/alignments/{tech}_{hap}_{name}_to_ref.bam", sample=wildcards.sample, hap=wildcards.hap, name=fill, tech=wildcards.tech)
	else:
		NAMES = glob_wildcards("results/{sample}/haplotype/{tech}_{hap}_{name}.list".format(tech=wildcards.tech, sample=wildcards.sample, hap=wildcards.hap, name='{name}')).name
		flist = expand("results/{sample}/alignments/{tech}_{hap}_{name}_to_ref.bam", sample=wildcards.sample, hap=wildcards.hap, name=NAMES,tech=wildcards.tech)
		return {'bams': flist,'flag': 'results/{sample}/haplotype/.{tech}_{hap}_split.done'}
		#return expand("results/{sample}/alignments/{tech}_{hap}_{name}_to_ref.bam", sample=wildcards.sample, hap=wildcards.hap, name=NAMES,tech=wildcards.tech)

#GET HAPLOTYPED METHYLATED BAM FILES
rule extract_haps:
	input:
		read_loc = 'results/{sample}/haplotype/{tech}_{hap}_{name}.list',
		flag = 'results/{sample}/haplotype/.{tech}_{hap}_split.done'
	output:
		just_names = 'results/{sample}/haplotype/{tech}_{hap}_{name}.names',
		bam = 'results/{sample}/methyl_bam/{tech}_{hap}_{name}.bam'
	resources:
		mem_mb = 12000,
		hrs = 12,
		threads = 1
	shell:
		'''
		awk -v OFS="\\t" '{{print $1}}' {input.read_loc} > {output.just_names}
		awk -v OFS="\\t" '{{print $2}}' {input.read_loc} | sort -u | xargs -i samtools view -b -o {output.bam} -N {output.just_names} {{}}
		'''
#GENERATE FASTQ FILES FROM BAM FILES WITH METHYLATED DATA
def getBAM(wildcards):
	if wildcards.hap == 'unphased':
		fofn = manifest_df.at[wildcards.sample,'bam_fofn']
		fofn_df = pd.read_csv(fofn, header=None, sep="\t", names=['file'])
		return fofn_df.at[int(wildcards.name),'file'] 
	else:
		return 'results/{sample}/methyl_bam/{tech}_{hap}_{name}.bam'
rule bam_to_fastq:
	input:
		bam = getBAM,
	output:
		fastq = 'results/{sample}/methyl_bam/{tech}_{hap}_{name}.fastq.gz'	
	resources:
		mem_mb = 8000,
		hrs = 12,
		threads = 8
    # added module load htslib/1.17-19-g07638e to get bgzip in PATH (DG, Nov 15, 2024).  conda's samtools didn't work
	shell:
		'''
		module load samtools/1.20 && module load htslib/1.17-19-g07638e1 && samtools fastq -T MM,ML -@{threads} {input.bam} | bgzip -c > {output.fastq}
		'''
rule copy_ref:
	input:
		ref =  lambda wildcards: manifest_df.at[wildcards.sample,'asm'],
	output:
		cp_ref = 'data/{sample}/align_to_asm.fa',
		cp_ref_fai = 'data/{sample}/align_to_asm.fa.fai'
	resources:
		mem_mb = 4000,
		hrs = 5,
		threads = 1
	shell:
		'''
		cp {input.ref} {output.cp_ref}
		samtools faidx {output.cp_ref}
		'''
#ALIGN SAMPLE METHYLATED FASTQ FILE TO REF
rule align_to_ref:
	input:
		ref = rules.copy_ref.output.cp_ref,
		fastq = 'results/{sample}/methyl_bam/{tech}_{hap}_{name}.fastq.gz'
	output:
		bam = 'results/{sample}/alignments/{tech}_{hap}_{name}_to_ref.bam'
	resources:
		mem_mb = 50000,
		hrs = 24,
		threads = 16
	params:
		sort_t = 2,
		tech = lambda wildcards: tech_dict[manifest_df.at[wildcards.sample,'tech_type']]
	# Without this, it doesn't work.  I think it allows the inherited snakemake_plus environment to be used
	# conda:
	# 	"envs/cpg.yaml"
	shell:
		'''
		meryl count k=15 output merylDB {input.ref} 
		mkdir -p rep_k15 && meryl print greater-than distinct=0.9998 merylDB > rep_k15/{wildcards.sample}_repetitive_k15.txt
		module load samtools/1.20 && winnowmap -W {wildcards.sample}_repetitive_k15.txt -y --eqx -ax {params.tech} -s 4000 -t 16 -I 10g {input.ref} {input.fastq} | samtools view -u -F 2308 - | samtools sort -o {output.bam} -
		'''
#COMBINE BAMS OF SAMPLE ALIGNED TO REF
rule combine_bams:
	input:
		unpack(find_reads),
	output:
		combined_bam = "results/{sample}/alignments/{tech}_{hap}_cpg.bam",
		combined_bai = "results/{sample}/alignments/{tech}_{hap}_cpg.bam.bai"
	resources:
		mem_mb = 8000,
		hrs = 24,
		threads = 12
	conda:
		"cpg2"
	shell:
		'''
		samtools merge -@{threads} -o {output.combined_bam} {input.bams}
		samtools index {output.combined_bam}
		'''
#CONVERT BAM FILE TO METHYL BED FILE

rule modbam2bed:
	input:
		ref = rules.copy_ref.output.cp_ref,
		bam = rules.combine_bams.output.combined_bam
	output:
		bed = 'results/{sample}/cpg/{sample}_{tech}_{hap}.bed'
	resources:
		mem_mb = 12000,
		hrs = 24,
		threads = 8
	conda:
		"cpg2"
	shell:
		'''
		modbam2bed -p {wildcards.sample}_{wildcards.tech}_{wildcards.hap} -m 5mC -t {threads} -e --aggregate --cpg {input.ref} {input.bam} > {output.bed}
		'''


#rule test:
#	input: 'test.txt'
#	output: touch('results/{sample}/cpg/{tech}_{hap}.bed')


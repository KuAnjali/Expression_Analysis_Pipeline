#!/usr/bin/python
'''
Function for tools used in the pipeline are defined here.
'''

import ConfigParser, os,sys
from os.path import dirname, realpath, sep, pardir
from functions import *

config_file  ='app.cfg'
config = ConfigParser.ConfigParser()
config.readfp(open(config_file))

def fastQC(outDir,ProjectID):
	commands=[]
	cmd1='python %s/1_genMat.py %s/%s/cfg %s %s'%(config.get('FASTQC_scripts','genMat'),outDir,ProjectID,40,outDir)
	commands.append(cmd1)
	cmd2='python %s/2_genFig.py %s/%s/cfg %s %s'%(config.get('FASTQC_scripts','genFig'),outDir,ProjectID,40,outDir)
	commands.append(cmd2)
	cmd3='python %s/3_genPdf.py %s/%s/cfg %s %s'%(config.get('FASTQC_scripts','genPdf'),outDir,ProjectID,40,outDir)
	commands.append(cmd3)
	cmd4='python %s/4_email.py %s/%s/cfg %s %s'%(config.get('FASTQC_scripts','email'),outDir,ProjectID,40,outDir)
	commands.append(cmd4)
	start_index=0
	while start_index < len(commands):
		sysVal_fastqc= osSys(str(commands[start_index]))
		start_index+=1		
        return  


def adapTrim(sampleDir,ClientID,R1,R2):
	adapt_cmd='java -Xmx8000M -jar %s/trimmomatic-0.36.jar PE %s %s %s/%s.R1_AT.PE.fq.gz %s/%s.R1_AT.SE.fq.gz %s/%s.R2_AT.PE.fq.gz %s/%s.R2_AT.SE.fq.gz -threads 32 ILLUMINACLIP:%s/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &>> %s/%s_trimmed_reads.log'%(config.get('tools','trimmomatic'),R1,R2,sampleDir,ClientID, sampleDir,ClientID,sampleDir,ClientID,sampleDir,ClientID, config.get('tools','trimmomatic'),sampleDir,ClientID)
	osSys(adapt_cmd)
	return

def contamination_removal(sampleDir,ClientID):
	contam_rm='%s/bowtie2 -x %s/hg19.contamination.fa -p 30 -N 1 -1 %s/%s.R1_AT.PE.fq.gz -2 %s/%s.R2_AT.PE.fq.gz --un-conc %s/%s.uncontaminated.fq  -S %s/%s.uncontaminated.sam &> %s/%s.contamination_removal.log'%(config.get('tools','bowtie2'),config.get('files','bowtie_hg19_index'),sampleDir,ClientID,sampleDir,ClientID,sampleDir,ClientID,sampleDir,ClientID,sampleDir,ClientID)
	osSys(contam_rm)
	return

def hisat2(sampleDir,ClientID):
	commands=[]
	cmd1='%s/hisat2 -p 32  --dta-cufflinks -x %s/genome -1 %s/%s.uncontaminated.1.fq -2 %s/%s.uncontaminated.2.fq -S %s/%s_HISAT2_aligned.sam --summary-file %s/%s_alignment_summary.txt'%(config.get('tools','hisat2'), config.get('files','hisat_hg19_index'),sampleDir,ClientID,sampleDir,ClientID,sampleDir,ClientID,sampleDir,ClientID)  
	commands.append(cmd1)
	cmd2='%s/samtools view -Sbuh -o %s/%s_HISAT2_aligned.bam %s/%s_HISAT2_aligned.sam'%(config.get('tools','samtools'),sampleDir,ClientID, sampleDir,ClientID)
	commands.append(cmd2)
	start_index=0
        while start_index < len(commands):
                hisat_cmd=osSys(str(commands[start_index]))
                start_index+=1
        return


def bam_sorting_indexing(sampleDir,ClientID):
	commands=[]
        sorting_cmd='%s/samtools sort -T temp -@ 5 -o %s/%s_HISAT2_aligned.sorted.bam -O bam %s/%s_HISAT2_aligned.bam'%(config.get('tools','samtools'),sampleDir,ClientID,sampleDir,ClientID)
        osSys(sorting_cmd)
        indexing_cmd='%s/samtools index %s/%s_HISAT2_aligned.sorted.bam'%(config.get('tools','samtools'),sampleDir,ClientID)
        osSys(indexing_cmd)
        return


def qualimap_read_align_dist(sampleDir,ClientID):
	cmd='%s/qualimap rnaseq -bam %s/%s_HISAT2_aligned.bam --java-mem-size=4G -gtf %s/NEW_Homo_sapiens.GRCh37.75.gtf -outdir %s/Qualimap/'%(config.get('tools','qualimap'),sampleDir, ClientID,config.get('files','hg19_gtf_star'),sampleDir)
	osSys(cmd)
	return

def rseqc_splice_junction(sampleDir, ClientID):
	cmd='python %s/junction_annotation.py -i %s/%s_HISAT2_aligned.bam -o %s/Splice_Junction/%s_junction_annot -r %s/NEW_Homo_sapiens.GRCh37.75.bed12 &> %s/Splice_Junction/%s_junction_annot_summary.txt'%(config.get('tools','rseqc'),sampleDir,ClientID,sampleDir,ClientID,config.get('files','human_hg19_bed'),sampleDir,ClientID)
	osSys(cmd)
	return


def cuffLinks(cufflink_path,sampleDir,ClientID):
	cufflink_cmd='%s/cufflinks -o %s  -p 15 --library-type fr-firststrand --max-bundle-length 6000000 --max-bundle-frags 10000000000000 -G %s/NEW_Homo_sapiens.GRCh37.75.gtf -L %s %s/%s_HISAT2_aligned.sorted.bam'%(config.get('tools','cufflink'),cufflink_path,config.get('files','hg19_gtf_star'), ClientID,sampleDir,ClientID)
	osSys(cufflink_cmd)
	return


def featureCount(sampleDir,ClientID):
	commands=[]
	feature_cmd='%s/featureCounts -Q 10 -s 0 -T 8 -a %s/NEW_Homo_sapiens.GRCh37.75.gtf -o %s/%s_FeatureCount.txt %s/%s_HISAT2_aligned.sorted.bam'%(config.get('tools','feature-count'),config.get('files','hg19_gtf_star'), sampleDir,ClientID, sampleDir,ClientID)
	commands.append(feature_cmd)
	extract_readcount="awk 'FNR > 2 {print $1, $NF}'  %s/%s_FeatureCount.txt >>  %s/%s_FeatureCount_final.txt"%(sampleDir,ClientID,sampleDir,ClientID)
	commands.append(extract_readcount)
	start_index=0
        while start_index < len(commands):
                sysVal_featureCount= osSys(str(commands[start_index]))
                start_index+=1
        return

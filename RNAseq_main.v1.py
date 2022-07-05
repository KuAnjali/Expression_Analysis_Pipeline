################################################
#                                              #
#                                              #
#       Author: Anjali Kumari                  #
#                                              #
#                                              #
###############################################

#!/usr/bin/python
import sys,re,os,csv,glob
import subprocess
from functions import *
import logging
import logging.config
from tools import *
from reportRmarkdown import *
import numpy as np
from stats_functions import *
import shutil


logging_file = 'logging.conf'
config_file  = 'app.cfg'

logging.config.fileConfig(logging_file)
logger=logging.getLogger(__name__)

#python RNA.py input.txt config.txt
if len(sys.argv)==3:
	master_config_dict={}
	master_config=open(sys.argv[2],"r")
	for line in master_config:
		line=line.strip("\n")
		line=line.split("\t")		
		master_config_dict[str(line[0])]=",".join(line[1:])
	mailID=str(master_config_dict["MAIL_ID"])

	input_file=open(sys.argv[1],'r').readlines()[1:]	
	if int(master_config_dict["FASTQC"])== 1:
		for line in input_file:
			line=line.strip("\n")
			line=line.split("\t")
			ProjectID=line[0]
			Medgenome_SampleID=line[1]
			ClientID=line[2]
			outDir=line[5]
			makedir='mkdir -p %s/%s'%(outDir,ProjectID)
			osSys(makedir)
			fastqc_path=outDir+'/'+ProjectID+'/FastQC'
			makedir='mkdir -p %s'%(fastqc_path)
			osSys(makedir)
			MAIL='MAIL\t'+mailID
			row='PE\t'+line[2]+'\t\t'+line[3]+'\t'+line[4]+'\n'
			headers=('PID\t'+ProjectID+'\nEID\tFastQC\n')
			cfg_file=outDir+'/'+ProjectID+'/cfg'
			cfg=open(outDir+'/'+ProjectID+'/cfg', 'a')
			cfg.write(row)	
			cfg.close()
		
		pre_postpend(cfg_file,headers,MAIL)
		logger.info("FastQC of the samples started!")
		fastQC(outDir,ProjectID)
		if str(os.path.exists(ProjectID+'_FastQC.pdf'))=='True':
			logger.info("FastQC of the samples completed and PDF has been sent to the provided mail ids!")
		else:
			logger.error("FastQC report has not been generated!")

	for line in input_file:
		line=line.strip("\n")
		line=line.split("\t")
		ProjectID=line[0]
		Medgenome_SampleID=line[1]
		ClientID=line[2]
		R1_FASTQ=line[3]
		R2_FASTQ=line[4]
		outDir=line[5]
		experi_group=line[6]
		makedir='mkdir -p %s/%s'%(outDir,ProjectID)
		osSys(makedir)
		sampleDir=outDir+'/'+ProjectID+'/'+ClientID
		makedir='mkdir -p %s'%(sampleDir)
		osSys(makedir)
		sampleInfo_file=outDir+'/'+ProjectID+'/sampleInfo.txt'
		sampleInfo=open(outDir+'/'+ProjectID+'/sampleInfo.txt','a')
		outLine=ClientID+'\t'+Medgenome_SampleID+'\t'+experi_group
		head=('Sample-Label\tMedgenomeID\tExperimental Group\n')
		sampleInfo.write(outLine+'\n')
		sampleInfo.close()
		if int(master_config_dict["FASTQC_SUMMARY"])==1:
			logger.info("FastQC summary of "+ClientID+" Started!")
			fastqc_stats(ProjectID,sampleDir,ClientID,outDir)
			logger.info("FastQC summary of "+ClientID+" Finished!")
						
		if int(master_config_dict["ADAP_TRIM"])==1:
			logger.info("Adapter Trimming of sample "+ClientID+" Started!")
			adapTrim(sampleDir,ClientID,R1_FASTQ,R2_FASTQ)
			logger.info("Adapter Trimming of sample "+ClientID+" Finished!")
		
		if int(master_config_dict["CONTAMINATION_REMOVAL"])==1:
			logger.info("Contamination Removal of sample "+ClientID+" Started!")
			contamination_removal(sampleDir,ClientID)
			logger.info("Contamination Removal of sample "+ClientID+" Finished!")
		
		if int(master_config_dict["ALIGNMENT"])==1:
			logger.info("Alignment of the sample "+ClientID+" Started!")
			hisat2(sampleDir,ClientID)	
			logger.info("Alignment of the sample "+ClientID+" Finished!")
			readCount_stats(outDir,ProjectID,sampleDir,ClientID)
			for f in os.listdir(sampleDir):
                		if f.endswith(('.sam','.fq.gz','.fq')):
					os.remove(os.path.join(sampleDir, f))

		if int(master_config_dict["BAM_PROCESSING"])==1:
			logger.info("Sorting and Indexing of BAM file "+ClientID+" Started!")
			bam_sorting_indexing(sampleDir,ClientID)
			logger.info("Finished Sorting and Indexing of BAM file "+ClientID+" !")

		if int(master_config_dict["SEQUENCE_ALIGNED_DISTRIBUTION"])==1:
			makedir='mkdir -p %s/Qualimap'%(sampleDir)
			osSys(makedir)
			logger.info("Aligned Read Distribution using Qualimap "+ClientID+" Started!")
			qualimap_read_align_dist(sampleDir,ClientID)
			logger.info("Aligned Read Distribution using Qualimap "+ClientID+" Finished!")
			qualiMap_stats(sampleDir,outDir,ProjectID,ClientID)
			shutil.rmtree(sampleDir+'/Qualimap')
		
		if int(master_config_dict["SPLICE_JUNCTION_DISTRIBUTION"])==1:
			makedir='mkdir -p %s/Splice_Junction'%(sampleDir)
			osSys(makedir)
			logger.info("Splice junction distribution using RSeQC "+ClientID+" Started!")
			rseqc_splice_junction(sampleDir, ClientID)
			logger.info("Splice junction distribution using RSeQC "+ClientID+" Finished!")	
			spliceJunction_stats(sampleDir,outDir,ProjectID,ClientID)
			shutil.rmtree(sampleDir+'/Splice_Junction')
	
	     	if int(master_config_dict["CUFFLINK"])==1:
			cufflink_path=sampleDir+'/CUFFLINK'
			makedir='mkdir -p %s'%(cufflink_path)
			osSys(makedir)
			logger.info("CUFFLINK of sample "+ClientID+" Started!")
			cuffLinks(cufflink_path,sampleDir,ClientID)
			logger.info("CUFFLINK of sample "+ClientID+" Finished!")
			
		if int(master_config_dict["FEATURE_COUNT"])==1:
			logger.info("FEATURE-COUNT of sample "+ClientID+" Started!")
			featureCount(sampleDir,ClientID)		
			logger.info("FEATURE-COUNT of sample "+ClientID+" Finished!")
			for f in os.listdir(sampleDir):
                        	if f.endswith(('FeatureCount.txt','FeatureCount.txt.summary')):
                        		os.remove(os.path.join(sampleDir, f))

			merge_file_input=open(outDir+'/'+ProjectID+'/merge_file_input.txt','a')
			merge_file_input.write(sampleDir+'/'+ClientID+'_FeatureCount_final.txt\t'+ClientID+'\t'+experi_group+'\n')
			merge_file_input.close()
		
	if int(master_config_dict["FEATURE_COUNT"])==1:
		if str(os.path.exists(outDir+'/'+ProjectID+'/merge_file_input.txt')) == 'True' and int(os.path.getsize(outDir+'/'+ProjectID+'/merge_file_input.txt')) > 0:
			deseq_input(outDir,ProjectID)
		
		if int(master_config_dict["DESEQ2"])==1:
			Rscript=open(outDir+'/'+ProjectID+'/Deseq2.R','w')
			Rscript.write('library("DESeq2")\n')
			merge_file_input=outDir+'/'+ProjectID+'/Deseq_input_merge_count.txt'
			Rscript.write('countData  = read.table("'+merge_file_input+'",header=TRUE,sep="\\t",row.names=1)\n')
			lst=[]
			merge_file_input=open(outDir+'/'+ProjectID+'/merge_file_input.txt','r').readlines()
			for line in merge_file_input:
				lst.append([ x for x in line.split()])
				group = [ x[2] for x in lst]
			condition = '","'.join(group)
			Rscript.write('condition <- c("'+condition+'")\n')
			Rscript.write("colData <- data.frame(row.names=colnames(countData), condition=factor(condition, levels=c('Control','Case')))\n")
                	Rscript.write('dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~condition)\n')
               		Rscript.write('dds <- DESeq(dds)\n')
                	Rscript.write('result <- results(dds, contrast=c("condition","Control", "Case"))\n')
                	Rscript.write('norm_counts <- counts(dds, normalized=TRUE)\n')
			Rscript.write('write.table(result,file="'+outDir+'/'+ProjectID+'/Deseq_result.txt"'+',sep="\\t")\n')
			Rscript.write('write.table(norm_counts,file="'+outDir+'/'+ProjectID+'/Deseq_normalised_result.txt"'+',sep="\\t", quote=FALSE)\n')
			Rscript.close()
			rcommand='Rscript %s/%s/Deseq2.R'%(outDir,ProjectID)
			osSys(rcommand)

		if int(master_config_dict["DESEQ2_GENE_FILTER"])==1:
			#Filtering Upregulated genes
			logger.info("Started filtering Upregulated genes using log2foldchange>=1 and padj<=0.05")
			upregulated_genes="awk  'BEGIN{OFS=\"\t\";FS=\"\t\"} $3>=1 && $7<=0.05 {print $0}' %s/%s/Deseq_result.txt | sed 's/\"//g'>> %s/%s/Upregulated_genes_deseq2.txt"%(outDir,ProjectID,outDir,ProjectID)
			subprocess.call(upregulated_genes, shell=True)
			if str(os.path.exists(outDir+'/'+ProjectID+'/Upregulated_genes_deseq2.txt')) == 'True' and int(os.path.getsize(outDir+'/'+ProjectID+'/Upregulated_genes_deseq2.txt')) > 0:
				logger.info("Filtered upregulated genes result is "+outDir+'/'+ProjectID+'/Upregulated_genes_deseq2.txt')
                                add_genename="awk 'FNR==NR{a[$2]=$1;next} ($1 in a) {print a[$1],$0}' /MGMSTAR2/SHARED/RESOURCES/RNAseq_pipeline_files/ENSGID_Genename.txt %s/%s/Upregulated_genes_deseq2.txt >%s/%s/Upregulated_genename_deseq2.txt"%(outDir,ProjectID,outDir,ProjectID)
                                subprocess.call(add_genename, shell=True)
				top50_up="sort -k 3nr %s/%s/Upregulated_genes_deseq2.txt  | head -50 > %s/%s/Top50_diffexp_Upregulated_genes.txt"%(outDir,ProjectID,outDir,ProjectID)
				subprocess.call(top50_up, shell=True)
				normalized_top50_up="awk 'FNR==NR{a[$1];next}($1 in a){print $0}' %s/%s/Top50_diffexp_Upregulated_genes.txt %s/%s/Deseq_normalised_result.txt > %s/%s/Top50_Upregulated_normalised_readcount.txt"%(outDir,ProjectID,outDir,ProjectID,outDir,ProjectID)
				subprocess.call(normalized_top50_up, shell=True)
				normalized_top50_up_genename="awk 'FNR==NR{a[$2]=$1;next} ($1 in a) {print a[$1],$0}' /MGMSTAR2/SHARED/RESOURCES/RNAseq_pipeline_files/ENSGID_Genename.txt %s/%s/Top50_Upregulated_normalised_readcount.txt > %s/%s/Top50_Upregulated_normalised_readcount_genename.txt"%(outDir,ProjectID,outDir,ProjectID)
				subprocess.call(normalized_top50_up_genename, shell=True)
			else:
				logger.info("No upregulated genes are found using padj parameter "+outDir+'/'+ProjectID+'/Upregulated_genes_deseq2.txt')
			
			#Filtering Downregulated genes
			logger.info("Started filtering Downregulated genes using log2foldchange<=-1 and padj<=0.05")
			downregulated_genes="awk  'BEGIN{OFS=\"\t\";FS=\"\t\"} $3<=-1 && $7<=0.05 {print $0}' %s/%s/Deseq_result.txt | sed 's/\"//g'>> %s/%s/Downregulated_genes_deseq2.txt"%(outDir,ProjectID,outDir,ProjectID)
                        subprocess.call(downregulated_genes, shell=True)			
			if str(os.path.exists(outDir+'/'+ProjectID+'/Downregulated_genes_deseq2.txt')) == 'True' and int(os.path.getsize(outDir+'/'+ProjectID+'/Downregulated_genes_deseq2.txt')) > 0:
                        	
			        logger.info("Filtered Downregulated genes result is "+outDir+'/'+ProjectID+'/Downregulated_genes_deseq2.txt')
				add_genename="awk 'FNR==NR{a[$2]=$1;next} ($1 in a) {print a[$1],$0}' /MGMSTAR2/SHARED/RESOURCES/RNAseq_pipeline_files/ENSGID_Genename.txt %s/%s/Downregulated_genes_deseq2.txt >%s/%s/Downregulated_genename_deseq2.txt"%(outDir,ProjectID,outDir,ProjectID)
				subprocess.call(add_genename, shell=True)
				top50_down="sort -nk3 %s/%s/Downregulated_genes_deseq2.txt  | head -50 > %s/%s/Top50_diffexp_Downregulated_genes.txt"%(outDir,ProjectID,outDir,ProjectID)
                                subprocess.call(top50_down, shell=True)
                                normalized_top50_up="awk 'FNR==NR{a[$1];next}($1 in a){print $0}' %s/%s/Top50_diffexp_Downregulated_genes.txt %s/%s/Deseq_normalised_result.txt > %s/%s/Top50_Downregulated_normalised_readcount.txt"%(outDir,ProjectID,outDir,ProjectID,outDir,ProjectID)
                                subprocess.call(normalized_top50_up, shell=True)
				normalized_top50_down_genename="awk 'FNR==NR{a[$2]=$1;next} ($1 in a) {print a[$1],$0}' /MGMSTAR2/SHARED/RESOURCES/RNAseq_pipeline_files/ENSGID_Genename.txt %s/%s/Top50_Downregulated_normalised_readcount.txt > %s/%s/Top50_Downregulated_normalised_readcount_genename.txt"%(outDir,ProjectID,outDir,ProjectID)
                                subprocess.call(normalized_top50_down_genename, shell=True)
                        else:
                                logger.info("No downregulated genes are found using padj parameter "+outDir+'/'+ProjectID+'/Downregulated_genes_deseq2.txt')
			
			top50_up_combine="cat %s/%s/Top50_Upregulated_normalised_readcount_genename.txt > %s/%s/Combined_top50_genes_normalised_readcount_genename.txt"%(outDir,ProjectID,outDir,ProjectID)
			subprocess.call(top50_up_combine, shell=True)
			top50_down_combine="grep -v '^GeneName' %s/%s/Top50_Downregulated_normalised_readcount_genename.txt >> %s/%s/Combined_top50_genes_normalised_readcount_genename.txt"%(outDir,ProjectID,outDir,ProjectID)
			subprocess.call(top50_down_combine, shell=True)
			head1="head -1 %s/%s/Deseq_normalised_result.txt"%(outDir,ProjectID)
			cmd_subprocess=subprocess.Popen(head1,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
	                cmd_out=cmd_subprocess.communicate()[0]
         	       	samplename=cmd_out.strip('\n').replace(' ','')
			firstline="GeneName\tENSGID\t"+samplename
			header(outDir+'/'+ProjectID+'/Combined_top50_genes_normalised_readcount_genename.txt',firstline)
			
			deseq_table(outDir,ProjectID)

		if int(master_config_dict["PATHWAY"])==1:
			#Pathway for Downregulated genes
			if str(os.path.exists(outDir+'/'+ProjectID+'/Downregulated_genename_deseq2.txt'))== 'True' and int(os.path.getsize(outDir+'/'+ProjectID+'/Downregulated_genes_deseq2.txt')) > 0:
				logger.info("Pathway analysis for Downregulated genes using Reactome_2016, WikiPathways_2016, KEGG_2016, GO_Molecular_Function_2018,GO_Cellular_Component_2018,GO_Biological_Process_2018, Pfam_InterPro_Domains Started!")
				pathway_down_R=open(outDir+'/'+ProjectID+'/Pathway_downregulated.R','w')
				pathway_down_R.write('library(enrichR)\n')
				pathway_down_R.write('dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","Reactome_2016","WikiPathways_2016", "KEGG_2016","Pfam_InterPro_Domains")\n')
				lst=[]
				downregulated_file=open(outDir+'/'+ProjectID+'/Downregulated_genename_deseq2.txt','r')
				for line in downregulated_file:
					line=line.strip("\n")
					gene=line.split(" ")[0]
					lst.append( gene)
				genes = '","'.join(lst)
				pathway_down_R.write('enriched <- enrichr(c("'+genes+'"), dbs)\n')
				pathway_down_R.write('GO_Biological <- enriched[["GO_Biological_Process_2018"]]\n')
				pathway_down_R.write('write.table(GO_Biological, file ="'+outDir+'/'+ProjectID+'/GO_Biological_Process_Downregulated.csv",sep="\\t")\n')
				pathway_down_R.write('GO_Molecular <- enriched[["GO_Molecular_Function_2018"]]\n')
				pathway_down_R.write('write.table(GO_Molecular, file ="'+outDir+'/'+ProjectID+ '/GO_Molecular_Function_Downregulated.csv",sep="\\t")\n')
				pathway_down_R.write('GO_Cellular <- enriched[["GO_Cellular_Component_2018"]]\n')
				pathway_down_R.write('write.table(GO_Cellular, file ="'+outDir+'/'+ProjectID+ '/GO_Cellular_Component_Downregulated.csv",sep="\\t")\n')
				pathway_down_R.write('Reactome <- enriched[["Reactome_2016"]]\n')
				pathway_down_R.write('write.table(Reactome, file = "'+outDir+'/'+ProjectID+ '/Reactome_Downregulated.csv",sep="\\t")\n')
				pathway_down_R.write('WikiPathways <- enriched[["WikiPathways_2016"]]\n')
				pathway_down_R.write('write.table(WikiPathways, file ="'+outDir+'/'+ProjectID+ '/WikiPathways_Downregulated.csv",sep="\\t")\n')
				pathway_down_R.write('Pfam <- enriched[["Pfam_InterPro_Domains"]]\n')
				pathway_down_R.write('write.table(Pfam, file = "'+outDir+'/'+ProjectID+ '/Pfam_InterPro_Domains_Downregulated.csv",sep="\\t")\n')
				pathway_down_R.write('KEGG <- enriched[["KEGG_2016"]]\n')
				pathway_down_R.write('write.table(KEGG, file = "'+outDir+'/'+ProjectID+ '/KEGG_Downregulated.csv",sep="\\t")\n')
				pathway_down_R.close()
				rcommand='Rscript %s/%s/Pathway_downregulated.R'%(outDir,ProjectID)
                        	osSys(rcommand)
			else:
				logger.info("No Downregulated genes are found with the cuttoff! Therefore No pathway analysis for Downregulated genes!")
			
			if str(os.path.exists(outDir+'/'+ProjectID+'/Upregulated_genename_deseq2.txt'))== 'True' and int(os.path.getsize(outDir+'/'+ProjectID+'/Upregulated_genes_deseq2.txt')) > 0:

                        	#Pathway for Upregulated genes
			        logger.info("Pathway analysis for Upregulated genes using Reactome_2016, WikiPathways_2016, KEGG_2016, GO_Molecular_Function_2018,GO_Cellular_Component_2018,GO_Biological_Process_2018, Pfam_InterPro_Domains Started!")
                                pathway_up_R=open(outDir+'/'+ProjectID+'/Pathway_upregulated.R','w')
                                pathway_up_R.write('library(enrichR)\n')
                                pathway_up_R.write('dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","Reactome_2016","WikiPathways_2016", "KEGG_2016","Pfam_InterPro_Domains")\n')
                                lst=[]
                                upregulated_file=open(outDir+'/'+ProjectID+'/Upregulated_genename_deseq2.txt','r')
                                for line in upregulated_file:
                                        line=line.strip("\n")
                                        gene=line.split(" ")[0]
                                        lst.append( gene)
                                genes = '","'.join(lst)
                                pathway_up_R.write('enriched <- enrichr(c("'+genes+'"), dbs)\n')
                                pathway_up_R.write('GO_Biological <- enriched[["GO_Biological_Process_2018"]]\n')
                                pathway_up_R.write('write.table(GO_Biological, file ="'+outDir+'/'+ProjectID+'/GO_Biological_Process_Upregulated.csv",sep="\\t")\n')
                                pathway_up_R.write('GO_Molecular <- enriched[["GO_Molecular_Function_2018"]]\n')
                                pathway_up_R.write('write.table(GO_Molecular, file ="'+outDir+'/'+ProjectID+ '/GO_Molecular_Function_Upregulated.csv",sep="\\t")\n')
                                pathway_up_R.write('GO_Cellular <- enriched[["GO_Cellular_Component_2018"]]\n')
                                pathway_up_R.write('write.table(GO_Cellular, file ="'+outDir+'/'+ProjectID+ '/GO_Cellular_Component_Upregulated.csv",sep="\\t")\n')
                                pathway_up_R.write('Reactome <- enriched[["Reactome_2016"]]\n')
                                pathway_up_R.write('write.table(Reactome, file = "'+outDir+'/'+ProjectID+ '/Reactome_Upregulated.csv",sep="\\t")\n')
                                pathway_up_R.write('WikiPathways <- enriched[["WikiPathways_2016"]]\n')
                                pathway_up_R.write('write.table(WikiPathways, file ="'+outDir+'/'+ProjectID+ '/WikiPathways_Upregulated.csv",sep="\\t")\n')
                                pathway_up_R.write('Pfam <- enriched[["Pfam_InterPro_Domains"]]\n')
                                pathway_up_R.write('write.table(Pfam, file = "'+outDir+'/'+ProjectID+ '/Pfam_InterPro_Domains_Upregulated.csv",sep="\\t")\n')
                                pathway_up_R.write('KEGG <- enriched[["KEGG_2016"]]\n')
                                pathway_up_R.write('write.table(KEGG, file = "'+outDir+'/'+ProjectID+ '/KEGG_Upregulated.csv",sep="\\t")\n')
                                pathway_up_R.close()
                                command='Rscript %s/%s/Pathway_upregulated.R'%(outDir,ProjectID)
                                osSys(command)
			else:
				logger.info("No Upregulated genes are found with the cuttoff! Therefore No pathway analysis for Upregulated genes!")
			go_combined_pathway_file(outDir, ProjectID)
			pathway_files(outDir, ProjectID)
	
	header(sampleInfo_file,head)	
	removal(outDir,ProjectID)
	if int(master_config_dict["REPORT"])== 1:
		logger.info("Started writing Rmarkdown Rmd file!")
		report(outDir,ProjectID)
		logger.info("Finished the analysis!")

else:
	logger.error("Please provide the correct number of inputs: python RNAseq_main.v1.py input.txt config.txt")	

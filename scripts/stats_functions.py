#!/usr/bin/python
import sys,re,os,csv,glob
import subprocess
import logging
import logging.config
from functions import *
import numpy as np
logging_file = 'logging.conf'
config_file  = 'app.cfg'

logging.config.fileConfig(logging_file)
logger=logging.getLogger(__name__)

def fastqc_stats(ProjectID,sampleDir,ClientID,outDir):
	files_list=['-LEFT-stats','-RIGHT-stats']
	dict={'Total Reads':[],'Average Quality': [], 'GC Percentage': [], 'Total Data Q30': [], 'Total Data Q20-Q30': [],'Total Data Q10-Q20':[], 'Average Read Length':[]}

	outFile=open(sampleDir+'/'+ClientID+'_FastQSummary.txt','w')
	j=1
	while len(files_list)>=j:
		for file in files_list:	
			R1=open(outDir+'/'+ProjectID+'/FastQC/'+ClientID+file,'r')
			for line in R1:
				line=line.strip("\n")
			
				#Total Number of Reads
				read_count_pat= re.search("Number of reads.+ \:(.+)",line)	
				if read_count_pat:
					dict.setdefault('Total Reads',[]).append(read_count_pat.group(1))
                                        j+=1

				#Average Base quality
				qual_pat= re.search("Average quality.+ \:(.+)",line)
				if qual_pat:
					dict.setdefault('Average Quality',[]).append(qual_pat.group(1))
					j+=1
				
				#Percent GC                  : 51.851953
				gc_pat = re.search("Percent GC.+ \:(.+)",line)
				if gc_pat:
					dict.setdefault('GC Percentage',[]).append(gc_pat.group(1))
                                        j+=1
				
				#% bases qual >= 30 
				Q30_pat = re.search("(% bases qual >= 30).+ \:(.+)",line)
				if Q30_pat:
					dict.setdefault('Total Data Q30',[]).append(Q30_pat.group(2))
					j+=1
				
				#% bases qual >= 20 and < 30
				Q20_Q30_pat = re.search("% bases qual >= 20 and < 30 \:(.+)",line)
				if Q20_Q30_pat:
					dict.setdefault('Total Data Q20-Q30',[]).append(Q20_Q30_pat.group(1))
                                        j+=1
				
				#% bases qual >= 10 and < 20
				Q10_Q20_pat = re.search("% bases qual >= 10 and < 20 \:(.+)",line)
				if Q10_Q20_pat:
					dict.setdefault('Total Data Q10-Q20',[]).append(Q10_Q20_pat.group(1))
                                        j+=1	
				#Average Read Length
				read_len_pat= re.search("Read length.+ \:(.+)",line)
				if read_len_pat:
					dict.setdefault('Average Read Length',[]).append(read_len_pat.group(1))
					j+=1
					
	summation('Total No. of Reads',dict['Total Reads'],outFile)
	average('Average Base Quality',dict['Average Quality'],outFile)
	average('GC Percentage', dict['GC Percentage'],outFile)
	average('Total Data >= Q30 %', dict['Total Data Q30'],outFile)
	average('Total Data >= Q20 and < 30 %', dict['Total Data Q20-Q30'],outFile)
	average('Total Data >= Q10 and < 20 %', dict['Total Data Q10-Q20'],outFile)
	average('Average Read Length',dict['Average Read Length'],outFile)
	outFile.close()
	fastq_summary_file=open(sampleDir+'/'+ClientID+'_FastQSummary.txt','r')
        rowList=[]
	rowList.append(ClientID)
        
        for line in fastq_summary_file:
        	line =line.strip("\n")
               	if line.startswith('Total No. of Reads'):
                	rowList.append(line.split(":")[1])
                if line.startswith('Average Base Quality'):
                	rowList.append(line.split(":")[1])
                if line.startswith('GC Percentage'):
                	rowList.append(line.split(":")[1])
                if line.startswith('Total Data >= Q30 %'):
                	rowList.append(line.split(":")[1])
                if line.startswith('Total Data >= Q20 and < 30 %'):
                	rowList.append(line.split(":")[1])
                if line.startswith('Total Data >= Q10 and < 20 %'):
                	rowList.append(line.split(":")[1])
                if line.startswith('Average Read Length'):
                	rowList.append(line.split(":")[1])
	headers=['Sample Name','Total No. of Reads','Average Base Quality','GC Percentage','Total Data >= Q30 %','Total Data >= Q20 and < 30 %','Total Data >= Q10 and < 20 %','Average Read Length']
	with open(outDir+'/'+ProjectID+'/FastQSummary.tmp','a') as csvfile:
        		writer = csv.writer(csvfile, delimiter='\t')
			writer.writerow(headers)
        		writer.writerow(rowList)
    	csvfile.close()
	
	lines_seen = set() # holds lines already seen
	outfile = open(outDir+'/'+ProjectID+'/FastQSummary.txt', "w")
	for line in open(outDir+'/'+ProjectID+'/FastQSummary.tmp', "r"):
    		if line not in lines_seen: # not a duplicate
        		outfile.write(line)
        		lines_seen.add(line)
	outfile.close()
	
def readCount_stats(outDir,ProjectID,sampleDir,ClientID):	
	outFile=open(sampleDir+'/'+ClientID+'_ReadCount_Summary.txt','w')
	#Read Count after adapter trimming #P27654_O51865/Control-K6-2/Control-K6-2_trimmed_reads.log
	alignPlot=open(sampleDir+'/'+ClientID+'_Alignment_Plot.txt','w')
	trim_Infile=open(sampleDir+'/'+ClientID+'_trimmed_reads.log','r')
	for line in trim_Infile:
		line=line.strip("\n")
		if line.startswith('Input'):
			#Input Read Pairs: 35619078 Both Surviving: 34369642 (96.49%)	
			pattern=re.search('Input Read Pairs: (.+) Both Surviving: (.+) \(.+\) Forward',line)
			Input_Reads=int(pattern.group(1))*2
			Reads_after_adap_trim=int(pattern.group(2))*2
			outFile.write("Total Reads :"+str(Input_Reads)+'\n')
			outFile.write("Read Count After Adapter Trimming :"+str(Reads_after_adap_trim)+'\n')
			adaptTrim_reads=int(Input_Reads)-int(Reads_after_adap_trim)
			alignPlot.write("Adapter Trimming :"+str(adaptTrim_reads)+'\n')
	#P27654_O51865/Control-K6-2/Control-K6-2.contamination_removal.log
				
	contam_Infile=open(sampleDir+'/'+ClientID+'.contamination_removal.log','r')
	for line in contam_Infile:
		line=line.strip("\n")
		#40645010 pairs aligned concordantly 0 times; of these:	
		pattern=re.search('(.+) pairs aligned concordantly 0 times; of these:',line)
		if pattern:
			Read_count_after_contam=int(pattern.group(1))*2
			outFile.write('Read Count After Contamination Removal :'+str(Read_count_after_contam)+'\n')
			#Contaminated reads = Reads after adapter trimming - Reads after contamination removal
			contaminated_Reads=int(int(Reads_after_adap_trim)-int(Read_count_after_contam))
			alignPlot.write("Contamination removal :"+str(contaminated_Reads)+'\n')
	#Control-K6-2_alignment_summary.txt
	align_Hisat2=open(sampleDir+'/'+ClientID+'_alignment_summary.txt','r')
	lines=[]
	for line in align_Hisat2:
		line=line.strip("\n")	
		pattern1=re.search('(.+) reads; of these:',line)
		if pattern1:
			lines.append(int(pattern1.group(1))*2)

		pattern2=re.search('(.+)% overall alignment rate',line)
		if pattern2:
			lines.append(float(pattern2.group(1)))
	aligned_reads=int((lines[0]*lines[1])/100)
	outFile.write("Reads Aligned :"+str(aligned_reads)+'\n')
	outFile.write("Alignment % :"+str(lines[1])+'\n')
	outFile.close()
	alignPlot.write("Reads Aligned :"+str(aligned_reads)+'\n')
	#Unligned Reads = Total Reads -(adaptTrim_reads+contaminated_Reads+aligned_reads)
	unaligned_Reads=int(Input_Reads)-int(int(adaptTrim_reads)+int(contaminated_Reads)+int(aligned_reads))
	alignPlot.write("Unaligned Reads :"+str(unaligned_Reads)+'\n')
	alignPlot.close()

	alignPlot_infile=open(sampleDir+'/'+ClientID+'_Alignment_Plot.txt','r')
	row_list=[]
        row_list.append(ClientID)
        for l in alignPlot_infile:
                l =l.strip("\n")
                if l.startswith('Adapter Trimming'):
                        row_list.append(l.split(":")[1])
                if l.startswith('Contamination removal'):
                        row_list.append(l.split(":")[1])
                if l.startswith('Reads Aligned'):
                        row_list.append(l.split(":")[1])
                if l.startswith('Unaligned Reads'):
                        row_list.append(l.split(":")[1])
        head=['SampleName','Adapter Reads','Contaminated Reads','Aligned Reads','Unaligned Reads']
        with open(outDir+'/'+ProjectID+'/Alignment_Plot.tmp','a') as alignfile:
                writer = csv.writer(alignfile, delimiter='\t')
                writer.writerow(head)
                writer.writerow(row_list)
        alignfile.close()

        l_seen = set() # holds lines already seen
        alignPlot_outFile = open(outDir+'/'+ProjectID+'/Alignment_Plot.txt', "w")
        for line in open(outDir+'/'+ProjectID+'/Alignment_Plot.tmp', "r"):
                if line not in l_seen: # not a duplicate
                        alignPlot_outFile.write(line)
                        l_seen.add(line)
        alignPlot_outFile.close()

	#Opening ReadCount summary file
	readcount_summary_file=open(sampleDir+'/'+ClientID+'_ReadCount_Summary.txt','r')
        rows_list=[]
        rows_list.append(ClientID)
        for line in readcount_summary_file:
        	line =line.strip("\n")
        	if line.startswith('Total Reads'):
     			rows_list.append(line.split(":")[1])
                if line.startswith('Read Count After Adapter Trimming'):
                	rows_list.append(line.split(":")[1])
                if line.startswith('Read Count After Contamination Removal'):
                	rows_list.append(line.split(":")[1])
                if line.startswith('Reads Aligned'):
                	rows_list.append(line.split(":")[1])
                if line.startswith('Alignment %'):
                	rows_list.append(line.split(":")[1]) 
	
	headers=['Sample Name','Total Reads','Read Count After Adapter Trimming','Read Count After Contamination Removal','Reads Aligned','Alignment']
 	with open(outDir+'/'+ProjectID+'/ReadCount_summary.tmp','a') as csvfile:
        	writer = csv.writer(csvfile, delimiter='\t')
		writer.writerow(headers)
        	writer.writerow(rows_list)
        csvfile.close()

	lines_seen = set() # holds lines already seen
        outfile = open(outDir+'/'+ProjectID+'/ReadCount_summary.txt', "w")
        for line in open(outDir+'/'+ProjectID+'/ReadCount_summary.tmp', "r"):
                if line not in lines_seen: # not a duplicate
                        outfile.write(line)
                        lines_seen.add(line)
        outfile.close()

	       
def qualiMap_stats(sampleDir,outDir,ProjectID,ClientID):
	seq_align_dist_file=open(sampleDir+'/Qualimap/rnaseq_qc_results.txt')
        output=[]
        output.append(ClientID)
        for line in seq_align_dist_file:
        	line=line.strip("\n")
                line=line.strip(" ")
                exon= re.search('exonic = (.+)\(.+\)',line)
                if exon:
                	output.append(exon.group(1))
		intron = re.search('intronic = (.+)\(.+\)',line)
                if intron:
                	output.append(intron.group(1))
		intergenic = re.search('intergenic = (.+)\(.+\)',line)
                if intergenic:
                	output.append(intergenic.group(1))
		overlapping = re.search('overlapping exon = (.+)\(.+\)',line)
                if overlapping:
                	output.append(overlapping.group(1))
	
	headers=['Sample Name','Exonic','Intronic','Intergenic','Overlapping Exon']
	with open(outDir+'/'+ProjectID+'/Sequence_alignment_dist.tmp','a') as csvfile:
        	writer = csv.writer(csvfile, delimiter="\t")
		writer.writerow(headers)
                writer.writerow(output)
        csvfile.close()

	lines_seen = set() # holds lines already seen
        outfile = open(outDir+'/'+ProjectID+'/Sequence_alignment_dist.txt', "w")
        for line in open(outDir+'/'+ProjectID+'/Sequence_alignment_dist.tmp', "r"):
                if line not in lines_seen: # not a duplicate
                        outfile.write(line)
                        lines_seen.add(line)
        outfile.close()
	
def spliceJunction_stats(sampleDir,outDir,ProjectID,ClientID):
	splice_junc_file=open(outDir+'/'+ProjectID+'/'+ClientID+'/Splice_Junction/'+ClientID+'_junction_annot_summary.txt','r')
	output=[]
	output.append(ClientID)
        for line in splice_junc_file:
                line=line.strip("\n")
                if line.startswith('Total splicing  Junctions:'):
                	line=line.replace("\t","")
                        output.append(line.split(':')[1])
                if line.startswith('Known Splicing Junctions:'):
                        line=line.replace("\t","")
                        output.append(line.split(':')[1])
                if line.startswith('Partial Novel Splicing Junctions:'):
                	line=line.replace("\t","")
                        output.append(line.split(':')[1])
                if line.startswith('Novel Splicing Junctions:'):
                	line=line.replace("\t","")
                        output.append(line.split(':')[1])

	headers=['SampleName','Total splicing Junctions','Known Splicing Junctions','Partial Novel Splicing Junctions','Novel Splicing Junctions']
	with open(outDir+'/'+ProjectID+'/Splice_Junction_Distribution.tmp','a') as csvfile:
        	writer = csv.writer(csvfile, delimiter="\t")
		writer.writerow(headers)
        	writer.writerow(output)
        csvfile.close()
	
	lines_seen = set() # holds lines already seen
        outfile = open(outDir+'/'+ProjectID+'/Splice_Junction_Distribution.txt', "w")
        for line in open(outDir+'/'+ProjectID+'/Splice_Junction_Distribution.tmp', "r"):
                if line not in lines_seen: # not a duplicate
                        outfile.write(line)
                        lines_seen.add(line)
        outfile.close()


def deseq_input(outDir,ProjectID):
        col_header = []
        genes = {}
	if str(os.path.exists(outDir+'/'+ProjectID+'/Deseq_input_merge_count.txt'))=='True':
        	os.remove(outDir+'/'+ProjectID+'/Deseq_input_merge_count.txt')

        outfile = open(outDir+'/'+ProjectID+'/Deseq_input_merge_count.txt','w')
        infile=open(outDir+'/'+ProjectID+'/merge_file_input.txt','r')
        for line in infile:
                line=line.strip(" ")
                filename,header,group= line.split('\t')
                data_f = open(filename,'r')
                col_header.append(header)
                #read file and add gene and counts to the dict
                for line in data_f:
                        gene,count = line.strip().split(' ')
                        if gene not in genes:
                                genes[gene] = [count]
                        else:
                                genes[gene].append(count)
                #important to close file
                data_f.close()
        infile.close()
        outfile.write('gene\t'+'\t'.join(col_header)+'\n')
        for gene in genes:
                data = genes[gene]
                #make sure each treatment has a count for this gene
                #this should catch most errors
                try:
                        assert len(data) == len(col_header)
                except AssertionError:
                        print "one of the treatment or genes is missing or extra"
                        print "data, found the problem here:"
                        print gene,data
                        print "while %s columns of treatments given" %len(col_header)
                        sys.exit()
                out_data = gene+'\t'+'\t'.join(data)+'\n'
                outfile.write(out_data)
        outfile.close()
	
def deseq_table(outDir,ProjectID):
	outFile=open(outDir+'/'+ProjectID+'/Differential_expressed_genes_summary.tmp','w')
	header="GeneName\tENSGID\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj"
	lst=[]
        merge_file_input=open(outDir+'/'+ProjectID+'/merge_file_input.txt','r').readlines()
        for line in merge_file_input:
		line=line.strip("\n")
		line=line.split("\t")
		lst.append(line[1])     
	outFile.write("Comparison :\t")
	outFile.write(",".join(str(item) for item in lst))

	if str(os.path.exists(outDir+'/'+ProjectID+'/Upregulated_genename_deseq2.txt')) == 'True' and int(os.path.getsize(outDir+'/'+ProjectID+'/Upregulated_genename_deseq2.txt')) > 0:
		totLines_up='cat %s/%s/Upregulated_genename_deseq2.txt | wc -l'%(outDir,ProjectID)
		cmd_subprocess=subprocess.Popen(totLines_up,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        	cmd_out=cmd_subprocess.communicate()[0]
        	num_upregulated_genes=cmd_out.strip('\n').replace(' ','')
		with open(outDir+'/'+ProjectID+'/Upregulated_genename_deseq2.txt','r+') as fp:
			content = fp.read()
			fp.seek(0, 0)
                	fp.write(header.rstrip('\r\n') + '\n' + content)
        	fp.close()
		outFile.write('\nNumber of Upregulated Genes :\t'+str(num_upregulated_genes)+'\n')

	else:
		num_upregulated_genes=0
		outFile.write('\nNumber of Upregulated Genes :\t'+str(num_upregulated_genes)+'\n')
	if str(os.path.exists(outDir+'/'+ProjectID+'/Downregulated_genename_deseq2.txt')) == 'True' and int(os.path.getsize(outDir+'/'+ProjectID+'/Downregulated_genename_deseq2.txt')) > 0:
                totLines_down='cat %s/%s/Downregulated_genename_deseq2.txt | wc -l'%(outDir,ProjectID)
		cmd_subprocess=subprocess.Popen(totLines_down,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                cmd_out=cmd_subprocess.communicate()[0]
                num_downregulated_genes=cmd_out.strip('\n').replace(' ','')
		with open(outDir+'/'+ProjectID+'/Downregulated_genename_deseq2.txt','r+') as fp:
			content = fp.read()
                        fp.seek(0, 0)
                        fp.write(header.rstrip('\r\n') + '\n' + content)
                fp.close()
		outFile.write('Number of Downregulated Genes :\t'+str(num_downregulated_genes)+'\n')

	else:
		num_downregulated_genes=0
		outFile.write('Number of Downregulated Genes :\t'+str(num_downregulated_genes)+'\n')
	outFile.close()
	
	finalOut=open(outDir+'/'+ProjectID+'/Differential_expressed_genes_summary.tmp','r')
	output=[]
        for line in finalOut:
                line=line.strip("\n")
                if line.startswith('Comparison :'):
                        line=line.replace("\t","")
                        output.append(line.split(':')[1])
                if line.startswith('Number of Upregulated Genes :'):
                        line=line.replace("\t","")
                        output.append(line.split(':')[1])
                if line.startswith('Number of Downregulated Genes :'):
                        line=line.replace("\t","")
                        output.append(line.split(':')[1])

	headers="Comparison\tNumber of Upregulated Genes\tNumber of Downregulated Genes"
	out=open(outDir+'/'+ProjectID+'/Differential_expressed_genes_summary.txt','w')
        out.write(headers+'\n')
	out.write("\t".join(str(item) for item in output))
	out.close()
		
def go_combined_pathway_file(outDir, ProjectID):
	GO_files=["GO_Biological_Process_Downregulated","GO_Cellular_Component_Downregulated","GO_Molecular_Function_Downregulated","GO_Biological_Process_Upregulated","GO_Cellular_Component_Upregulated","GO_Molecular_Function_Upregulated"]
	for filename in GO_files:
		inFile=outDir+'/'+ProjectID+'/'+filename+'.csv'
		if str(os.path.exists(inFile)) == 'True' and int(os.path.getsize(inFile)) > 0:	
			list1=str(filename).split("_")
			gene_set_library=list1[0]+" "+list1[1]+" "+list1[2]
			geneset=list1[3]
			go_combined_pathway_file=open(outDir+'/'+ProjectID+'/go_combined_pathway_'+geneset+'.tmp','a')
			for line in open(inFile,'r').readlines()[1:]:
				line=line.strip("\n")
				line=line.replace("\"","")
				line=line.split("\t")
				#regulation of endothelial cell chemotaxis to fibroblast growth factor (GO:2000544)	-3.73592390156943	11.1577130254489	0.713448947557876	0.0504586685185959	1/8	['FGF16']	GO Biological Process	Downregulated
				firstline= "Term"+'\t'+"Z.score"+'\t'+"Combined.Score"+'\t'+"Adjusted.P.value"+'\t'+"P.value"+'\t'+"Overlapping"+'\t'+"Genes"+'\t'+"gene_set_library"+'\t'+"geneset" 
				Z_Score=round(float(line[7]),2)
				Combined_score=round(float(line[8]),2)
				Padj=round(float(line[4]),3)
				Pvalue=round(float(line[3]),3)
				genes=line[9].split(";")
				go_combined_pathway_file.write(line[1]+'\t'+str(Z_Score)+'\t'+str(Combined_score)+'\t'+str(Padj)+'\t'+str(Pvalue)+'\t'+line[2]+'\t'+str(genes)+'\t'+str(gene_set_library)+'\t'+str(geneset)+'\n')
			
			go_combined_pathway_file.close()
			
			with open(outDir+'/'+ProjectID+'/go_combined_pathway_'+geneset+'.tmp', 'r+') as f:
                		content = f.read()
                		f.seek(0, 0)
                		f.write(firstline.rstrip('\r\n') + '\n' + content)
			f.close()
			
		else:
			continue

	lines_seen_up = set() # holds lines already seen
        outUp = open(outDir+'/'+ProjectID+'/go_combined_pathway_Upregulated.tsv', "w")
        for line in open(outDir+'/'+ProjectID+'/go_combined_pathway_Upregulated.tmp', "r"):
                if line not in lines_seen_up: # not a duplicate
                        outUp.write(line)
                        lines_seen_up.add(line)
        outUp.close()

	lines_seen_down = set() # holds lines already seen
        outDown = open(outDir+'/'+ProjectID+'/go_combined_pathway_Downregulated.tsv', "w")
        for line in open(outDir+'/'+ProjectID+'/go_combined_pathway_Downregulated.tmp', "r"):
                if line not in lines_seen_down: # not a duplicate
                        outDown.write(line)
                        lines_seen_down.add(line)
        outDown.close()
	return

def pathway_files(outDir, ProjectID):
	FILES=['GO_Biological_Process_Downregulated','GO_Biological_Process_Upregulated','GO_Molecular_Function_Downregulated','GO_Molecular_Function_Upregulated','GO_Cellular_Component_Downregulated','GO_Cellular_Component_Upregulated','KEGG_Downregulated','KEGG_Upregulated','Reactome_Downregulated','Reactome_Upregulated','WikiPathways_Downregulated','WikiPathways_Upregulated','Pfam_InterPro_Domains_Upregulated','Pfam_InterPro_Domains_Downregulated']
	for filename in FILES:
		inFile=outDir+'/'+ProjectID+'/'+filename+'.csv'
		if str(os.path.exists(inFile)) == 'True' and int(os.path.getsize(inFile)) > 0:
			times=filename.count("_")
			
			if times==1:
                        	list1=str(filename).split("_")
                        	gene_set_library=list1[0]
                        	geneset=list1[1]
				pathway_file=open(outDir+'/'+ProjectID+'/'+gene_set_library+'_pathway_'+geneset+'.tmp','a')
                        	for line in open(inFile,'r').readlines()[1:]:
                                	line=line.strip("\n")
                                	line=line.replace("\"","")
                                	line=line.split("\t")
                                	firstline= "Term"+'\t'+"Z.score"+'\t'+"Combined.Score"+'\t'+"Adjusted.P.value"+'\t'+"P.value"+'\t'+"Overlapping"+'\t'+"Genes"+'\t'+"gene_set_library"+'\t'+"geneset"
                                	Z_Score=round(float(line[7]),2)
                                	Combined_score=round(float(line[8]),2)
                                	Padj=round(float(line[4]),3)
                                	Pvalue=round(float(line[3]),3)
					genes=line[9].split(";")
                                	pathway_file.write(line[1]+'\t'+str(Z_Score)+'\t'+str(Combined_score)+'\t'+str(Padj)+'\t'+str(Pvalue)+'\t'+line[2]+'\t'+str(genes)+'\t'+str(gene_set_library)+'\t'+str(geneset)+'\n')

                        	pathway_file.close()

                        	with open(outDir+'/'+ProjectID+'/'+gene_set_library+'_pathway_'+geneset+'.tmp', 'r+') as f:
                        		content = f.read()
                                	f.seek(0, 0)
                                	f.write(firstline.rstrip('\r\n') + '\n' + content)
                        	f.close()
				lines_seen = set() # holds lines already seen
        			out = open(outDir+'/'+ProjectID+'/'+gene_set_library+'_pathway_'+geneset+'.tsv', "w")
        			for line in open(outDir+'/'+ProjectID+'/'+gene_set_library+'_pathway_'+geneset+'.tmp', "r"):
                			if line not in lines_seen: # not a duplicate
                        			out.write(line)
                        			lines_seen.add(line)
        			out.close()


			elif times==3:
				#print filename
				list1=str(filename).split("_")
				gene_set_library=list1[0]+'_'+list1[1]+'_'+list1[2]
				geneset=list1[3]
                        	pathway_file=open(outDir+'/'+ProjectID+'/'+gene_set_library+'_'+geneset+'.tmp','a')
                        	for line in open(inFile,'r').readlines()[1:]:
                                	line=line.strip("\n")
                                	line=line.replace("\"","")
                                	line=line.split("\t")
                                	firstline= "Term"+'\t'+"Z.score"+'\t'+"Combined.Score"+'\t'+"Adjusted.P.value"+'\t'+"P.value"+'\t'+"Overlapping"+'\t'+"Genes"+'\t'+"gene_set_library"+'\t'+"geneset"
                                	Z_Score=round(float(line[7]),2)
                                        Combined_score=round(float(line[8]),2)
                                        Padj=round(float(line[4]),3)
                                        Pvalue=round(float(line[3]),3)
					genes=line[9].split(";")
	                               	pathway_file.write(line[1]+'\t'+str(Z_Score)+'\t'+str(Combined_score)+'\t'+str(Padj)+'\t'+str(Pvalue)+'\t'+line[2]+'\t'+str(genes)+'\t'+str(gene_set_library)+'\t'+str(geneset)+'\n')

                        	pathway_file.close()

                        	with open(outDir+'/'+ProjectID+'/'+gene_set_library+'_'+geneset+'.tmp', 'r+') as f:
                        		content = f.read()
                                	f.seek(0, 0)
                                	f.write(firstline.rstrip('\r\n') + '\n' + content)
                        	f.close()
				lines_seen = set() # holds lines already seen
        			out = open(outDir+'/'+ProjectID+'/'+gene_set_library+'_'+geneset+'.tsv', "w")
        			for line in open(outDir+'/'+ProjectID+'/'+gene_set_library+'_'+geneset+'.tmp', "r"):
                			if line not in lines_seen: # not a duplicate
                        			out.write(line)
                        			lines_seen.add(line)
        			out.close()

		else:
			continue

	return
		

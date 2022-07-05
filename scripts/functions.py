#!/usr/bin/python
'''
Functions for small task are defined here! 
'''

import os, sys, glob
import subprocess
from datetime import datetime
import numpy as np
import shutil

def fileExists(R1,R2):
        flag=''
        file_exists_R1= os.path.exists(R1)
        file_exists_R2= os.path.exists(R2)
        if str(file_exists_R1)!='True' or str(file_exists_R2)!='True' :
                flag='False'
        else:
                flag='True'
        return flag


def osSys(cmd):
        flag=''
        sysVal=os.system(cmd)
        
        if sysVal==0:
                flag='True'
        else:
                flag='False'
        return flag

def pre_postpend(filename,headers,MAIL):
        with open(filename, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                f.write(headers.rstrip('\r\n') + '\n' + content)
        f.close()
        with open(filename,'a') as f:
                f.write(MAIL)
        f.close()

def header(filename,value):
	with open(filename, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                f.write(value.rstrip('\r\n') + '\n' + content)
        f.close()


def average(string,list,outfile):
	list= np.array(list).astype(np.float)
	list = np.mean(list)
	outfile.write('\n'+string+' : '+str(list))
	return

def summation(string,list,outfile):
	list= np.array(list).astype(np.int)
        list = np.sum(list,dtype=np.int32)
        outfile.write('\n'+string+' : '+str(list))
        return

def removal(outDir,ProjectID):
	for f in os.listdir(outDir+'/'+ProjectID):
		if f.endswith(('.csv','.tmp')):
			os.remove(os.path.join(outDir+'/'+ProjectID, f))
			
	for f in os.listdir(outDir+'/'+ProjectID+'/FastQC'):
		if f.endswith(('-bd','-gc','-qd','.png','-bd-pdf','-gc.pdf','-qd.pdf','.aux','.tex','.toc','.log','.out')):
			os.remove(os.path.join(outDir+'/'+ProjectID+'/FastQC', f))
	
	return


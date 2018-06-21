#!/usr/bin/env python3

################################################################
# Program: ALSgeneScanner variant score-DNAscan
# Version 0.1
# Author: Alfredo Iacoangeli (alfredo.iacoangeli@kcl.ac.uk)
#################################################################

import sys

variant_file = sys.argv[1]
out_file =  sys.argv[2]


w=open('%s' %(out_file), 'w')

a=open('%s' %(variant_file))

a_lines=a.readlines()

w.write(a_lines[0].split('\t')[0]+"\t"+a_lines[0].split('\t')[1]+"\t"+a_lines[0].split('\t')[2]+"\t"+a_lines[0].split('\t')[3]+"\t"+a_lines[0].split('\t')[4]+"\t"+a_lines[0].split('\t')[5]+"\t"+a_lines[0].split('\t')[6]+"\t"+a_lines[0].split('\t')[7]+"\t"+a_lines[0].split('\t')[8]+"\tTotal_score\t"+a_lines[0].split('\t')[-6]+"\t"+a_lines[0].split('\t')[-5]+"\t"+"\t"+a_lines[0].split('\t')[11]+"\t"+a_lines[0].split('\t')[13]+"\t"+a_lines[0].split('\t')[15]+"\t"+a_lines[0].split('\t')[17]+"\t"+a_lines[0].split('\t')[19]+"\t"+a_lines[0].split('\t')[21]+"\t"+a_lines[0].split('\t')[23]+"\t"+a_lines[0].split('\t')[25]+"\t"+a_lines[0].split('\t')[31]+"\t"+a_lines[0].split('\t')[33]+"\t"+a_lines[0].split('\t')[28]+"\n")

a_lines.pop(0)       

for i in a_lines:
     a=0
     b=0
     c=0
     d=0
     e=0
     f=0
     g=0
     h=0
     j=0
     l=0
     o=0
     
     
     counter=0
     if i.split('\t')[11]=="D":
              counter+=1
              a=1
     if i.split('\t')[13]=="D":
              counter+=1
              b=1
     if i.split('\t')[13]=="P":
              counter+=1
              b=1
     if i.split('\t')[15]=="D":
              counter+=1
              c=1
     if i.split('\t')[15]=="P":
              counter+=1
              c=1
     if i.split('\t')[17]=="D":
              counter+=1
              d=1
     if i.split('\t')[19]=="A":
              counter+=1
              e=1
     if i.split('\t')[19]=="D":
              counter+=1
              e=1
     if i.split('\t')[21]=="H":
              counter+=1
              f=1
     if i.split('\t')[21]=="M":
              counter+=1
              f=1
     if i.split('\t')[23]=="D":
              counter+=1
              g=1		  
     if i.split('\t')[25]=="D":
              counter+=1
              h=1
     if i.split('\t')[31]=="D":
              counter+=1
              j=1
     if i.split('\t')[33]=="D":
              counter+=1
              l=1
     if i.split('\t')[28] != '.' and i.split('\t')[28] != '' :
	 
     	if float(i.split('\t')[28])>=15:
              counter+=1
              o=1
     w.write(i.split('\t')[0]+"\t"+i.split('\t')[1]+"\t"+i.split('\t')[2]+"\t"+i.split('\t')[3]+"\t"+i.split('\t')[4]+"\t"+i.split('\t')[5]+"\t"+i.split('\t')[6]+"\t"+i.split('\t')[7]+"\t"+i.split('\t')[8]+"\t"+str(counter)+'\t'+i.split('\t')[-18]+"\t"+i.split('\t')[-17]+"\t"+str(a)+"\t"+str(b)+"\t"+str(c)+"\t"+str(d)+"\t"+str(e)+"\t"+str(f)+"\t"+str(g)+"\t"+str(h)+"\t"+str(j)+"\t"+str(l)+"\t"+str(o)+"\n")
	 
w.close()





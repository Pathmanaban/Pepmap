
"""
Pepmap.py

Created by Pathmanaban Ramasamy on 18th Sep 2018
Copyright (c) 2018 Pathmanaban Ramasamy. All rights reserved.

"""

import csv
import os
import difflib
from Bio import SeqIO
from collections import defaultdict


def UPseq():
	records = list(SeqIO.parse("uniprot_sprot.fasta", "fasta"))
	return records
UPseqlis=UPseq()				
	
def UPseqfetch(UPseqlis,accsn,pep):

	for record in UPseqlis:
			recid=record.id.split('|')[1]
			if str(recid)==str(accsn).strip():
				
				upseq=record.seq
				
				newpos1= [n for n in xrange(len(upseq)) if upseq.find(pep, n) == n]
				if newpos1:
					newpos1=[x+1 for x in newpos1]
					
					return ("No",newpos1)
				
				
def extract_mod():
	with open('unimodptms.txt','r') as infilex:
		next(infilex,None)
		modific=[]
		for line in infilex:
			line=line.split(',')
			mod=line[0].split('=')[1].lower()
			if mod not in modific:
				modific.append(mod)
			
		
			
		return modific
		
	

modifi=extract_mod()

def map_and_locate():
	
	infile=raw_input("Enter the file name to process: ")
	outfile=raw_input("Enter the outfile name: ")
	propepindx,pepindx,accession = input("Enter column numbers of peptide, modified peptide and protein accession sepereated by comma: ")
	with open(infile,'r') as infile1,open(outfile,'w') as outfile:
			writer1 = csv.writer(outfile,delimiter='\t')
			
			next(infile1,None)
			decoy=['_crap','Random_']
			an,bn=0,0
			for line1 in infile1:
				an+=1
				
				line=line1.split('\t')
				protpep=line[propepindx]
				pep=line[pepindx]
				
				if '_HUMAN' in line[accession] and '_crap' not in line[accession] and 'Random_' not in line[accession]:
					
					if '|' in line[accession]:
						bn+=1
						accsn=line[accession].split('|')[1].strip()
						
						
						pep_start1=UPseqfetch(UPseqlis,accsn,protpep)
						if pep_start1:
							pep_start=pep_start1[1]
							mismatch=pep_start1[0]
					
							
							if pep_start:
								
								for ev_p_start in pep_start:
									ev_p_start=int(ev_p_start)
									peplis=[]
									modlis=[]
						
									if len(line[propepindx]) != len(line[pepindx]):
											
											for s in modifi:
											
												if s in pep.lower() and s not in modlis:
													
													modlis.append(s)
											
											for mod in modlis:		
												mod_pos_pep = [n for n in xrange(len(pep)) if pep.lower().find(mod, n) == n]
												
												for x in mod_pos_pep:
													peplis.append(mod+'__'+str(x+ev_p_start))
												
												peplis=sorted(peplis, key=lambda x: int(x.split('__')[1]))
												
											if len(peplis) >1:
												peplis1=[x.split('__')[0] for x in peplis]
												newpos=[]
												lis2=[len(x.split('__')[0]) for x in peplis]
												lis1=[int(x.split('__')[1]) for x in peplis]
												for elem in lis1:
													elem1=elem-sum(lis2[:lis1.index(elem)])
													newpos.append(elem1)
												
												peplis2=[]
												for a,b in enumerate(peplis1):
													
													elem2=b+'__'+str(newpos[a])
													peplis2.append(elem2)
												
												writer1.writerow([line1,ev_p_start,','.join(peplis2),mismatch])
											elif len(peplis)==1:
												writer1.writerow([line1,ev_p_start,','.join(peplis),mismatch])
									else:
										writer1.writerow([line1,ev_p_start,"Unmodified",mismatch])
	print "\n"									
	print "Results are wirtten to: ", outfile
	print "\n"
	print "*****************************************************"
	print "              Thanks for using Pepmap                "
	print "*****************************************************"

if __name__ == "__main__":
	print "\n"
	print "************************************************"
	print "              Pepmap version 1.0                "
	print "************************************************"
	print "\n"
	map_and_locate()				
	



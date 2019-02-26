#!/usr/bin/python
#usage: python script.py file.fasta
import sys
import re 
filefasta = sys.argv[1]
 
try:
	f = open(filefasta)
except IOError:
	print ("File doesn't exist!")
 


def translate(seq): 
	   
	table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				  
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	} 
	protein ="" 
	if len(seq)%3 == 0: 
		for i in range(0, len(seq), 3): 
			codon = seq[i:i + 3] 
			protein+= table[codon] 
	protein = protein[:-1]
	return protein 


 

def inversComplement(input):
	output = ''
	for letter in input:
		letter = letter.upper()

		if letter == 'A':
			output += 'T'
		elif letter == 'T':
			output += 'A'
		elif letter == 'G':
			output += 'C'
		elif letter == 'C':
			output += 'G'
		elif letter == 'N':
			output += 'N'
		else:
			output += '-'
			print letter
	return(output[::-1])
 
	



# print seq
# print seqAnti

# print seq
def getORF(seq):
	dicPath = dict()
	for i in range (0, len(seq) - 3):
		tmpSeq = seq[i:]
		p = re.finditer(r'(ATG((?:.{3})+?)T(?:AG|AA|GA))', tmpSeq)
		if p:
			for e in p:
				if not 'N' in e.group() and len(e.group()) >= 100:
					aa = translate(e.group())
					if not '_' in aa:
						pos = 1 + seq.index(e.group())
						frame = pos % 3
						if frame == 0:
							frame = 3			
						dicPath[e.group()] = str(pos) + "\t" + str(len(e.group())) + "\t" + str(frame)
			p = None


	return dicPath


# print "############################"
# 
# 
# dicPath = dict()
# for i in range (0, len(seqAnti) - 3):
#	 tmpSeq = seqAnti[i:]
#	 p = re.finditer(r'(ATG((?:.{3})+?)T(?:AG|AA|GA))', tmpSeq)
#	 if p:
#		 for e in p:
#			 pos = 1 + len(seq)-(seqAnti.index(e.group()) + len(e.group()))
#			 frame = pos % 3 
#			 if frame == 0:
#				 frame = 3
# 
#			 dicPath[e.group()] = str(pos) + "\t-\t" + str(len(e.group())) + "\t" + str(frame)
#		 p = None
#		 
# 
# for x, y in dicPath.items():
#   print(x, y) 



# for e in p:
#	 print e
	# print e.start()
#	 print len(e[0])
#	 print e[0]
# 
# print seqAnti
# print 'anti'
# p = re.findall(r'(ATG((?:.{3})+?)T(?:AG|AA|GA))', seqAnti)
# for e in p:
#	 
#	 print len(e[0])
#	 print e[0]
# 
# (?=A[TU]G((?:.{3})+?)[TU](?:AG|AA|GA))
# print p
# m = re.compile("^(?:gau|gac)(?:auu|auc|aua)(?:ggu|ggc|gga|ggg)(?:ggu|gac)$")





# 
# 
# 
# 
# 
# def orfFINDER(dna,frame):
#	  
#	 stop_codons = ['tga', 'tag', 'taa']
#	 start_codon = ['atg']
#	 start_positions = []
#	 stop_positions = []
#	 num_starts=0
#	 num_stops=0
#	  
#	  
#	 for i in range(frame,len(dna),3):
#		 codon=dna[i:i+3].lower()
#		 if codon in start_codon:
#			 start_positions += str(i+1).splitlines()
#		 if codon in stop_codons:
#			 stop_positions += str(i+1).splitlines()
#	  
#	 for line in stop_positions:
#		 num_stops += 1
#	  
#	 for line in start_positions:
#		 num_starts += 1
#	  
#	 orffound = {}
#	  
#	 if num_stops >=1 and num_starts >=1: #first statment: the number of stop codons and start condos are greater than or equal to 1;
#		  
#		 orfs = True
#		 stop_before = 0
#		 start_before = 0
#		  
#		 if num_starts > num_stops:
#			 num_runs = num_starts
#		 if num_stops > num_starts:
#			 num_runs = num_stops
#		 if num_starts == num_stops:
#			 num_runs = num_starts
#			  
#		 position_stop_previous = 0
#		 position_start_previous = 0
#		 counter = 0
#		  
#		 for position_stop in stop_positions:
#			  
#			 position_stop = int(position_stop.rstrip()) + 2
#			  
#			 for position_start in start_positions:
#				  
#				 position_start = position_start.rstrip()
#				  
#				 if int(position_start) < int(position_stop) and int(position_stop) > int(position_stop_previous) and int(position_start) > int(position_stop_previous):
#					  
#					 counter += 1
#					 nameorf = "orf"+str(counter)
#					 position_stop_previous += int(position_stop) - int(position_stop_previous)
#					 position_start_previous += int(position_start) - int(position_start_previous)
#					 sizeorf = int(position_stop) - int(position_start) + 1
#					  
#					 orffound[nameorf] = position_start,position_stop,sizeorf,frame
#						  
#				 else:
#					  
#					 pass
#	  
#	 else:
#		  
#		 orfs = False
#  
#	 return orffound
# #FUNCTION END
#  
#  
# #READ FASTA FILE AND SAVE HEADERS AND SEQUENCES IN A DICTIONARY
# seqs={}
#  

def writeToFile(head, dicORF, ori):
	for x, y in dicORF.items():
		print head + "\t" + ori + "\t" + x + "\t" + y + "\t" + translate(x)


seqs=dict()
for line in f:
	line = line.rstrip()
	if line[0] == '>':
		words=line.split() 
		name=words[0][1:]
		seqs[name]=''
	else:
		seqs[name] = seqs[name] + line
		# print seqs[name]
#DEFINE FRAME TO FIND ORF
#if frame = 0, start from the first position in the sequence
frame=0
 
#EXECUTE THE ORFFINDER FUNCTION
print "ID\tstrand\tORF\tstart_pos\tORF_length\tframe\tAA"
for i in seqs.items():
	header= i[0]
	seq = i[1]
	seqAnti=inversComplement(seq)
	
	dicFrw = getORF(seq)
	writeToFile(header, dicFrw, "+")
	dicRvs = getORF(seqAnti)
	writeToFile(header, dicRvs,  "-")
 

f.close()


# 

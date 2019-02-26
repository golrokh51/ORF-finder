#!/usr/bin/python
#usage: python script.py file.fasta
import sys
import re 
filefasta = sys.argv[1]
ntSize = 3
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
	if len(seq) % ntSize == 0: 
		for i in range(0, len(seq), ntSize): 
			codon = seq[i:i + ntSize] 
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
	for i in range (0, len(seq) - ntSize):
		tmpSeq = seq[i:]
		p = re.finditer(r'(ATG((?:.{3})+?)T(?:AG|AA|GA))', tmpSeq)
		if p:
			for e in p:
				if not 'N' in e.group() and len(e.group()) >= 100:
					aa = translate(e.group())
					if not '_' in aa:
						pos = 1 + seq.index(e.group())
						frame = pos % ntSize
						if frame == ntSize:
							frame = ntSize			
						dicPath[e.group()] = str(pos) + "\t" + str(len(e.group())) + "\t" + str(frame)
			p = None


	return dicPath




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
print "#ID\tstrand\tORF\tstart_pos\tORF_length\tframe\tAA"
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

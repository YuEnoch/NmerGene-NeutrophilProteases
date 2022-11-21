'''Example of a python module.
'''

#This is the one letter universal translation table
#It handles cases of DNA ambiguity where the encoded amino acid is unambiguous.
#You need to deal with the missing cases where ambiguity codes would result in
#an ambiguous amino acid assignment. It is suggested that you use 'X' in these
#cases as this is the standard character for an unknown amino acid.
#Only Y (pyrimidine), R (purine) and N (any) degeneracy symbols are handled at
#this time. (need to add M,K,W,S,B,D,H,V where appropirate)
#Stop codons are symbolized as X
#Reassign TAA, TAG, TAR and TGA to change the stop codon sybmol if desired.

transTab1L = {
'TTT': 'F', 'TTC': 'F', 'TTY': 'F', 'TTA': 'L', 'TTG': 'L', 'TTR': 'L', 
'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TCN': 'S', 'TCY': 'S', 'TCR': 'S', 
'TAT': 'Y', 'TAC': 'Y', 'TAY': 'Y', 'TAA': 'X', 'TAG': 'X', 'TAR': 'X', 
'TGT': 'C', 'TGC': 'C', 'TGY': 'C', 'TGA': 'X', 'TGG': 'W', 
'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CTY': 'L', 'CTR': 'L', 'CTN': 'L',
	'YTG': 'L', 'YTA': 'L', 
'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCY': 'P', 'CCR': 'P', 'CCN': 'P', 
'CAT': 'H', 'CAC': 'H', 'CAY': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CAR': 'Q', 
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'CGY': 'R', 'CGR': 'R', 'CGN': 'R', 
'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATY': 'I', 'ATG': 'M', 
'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ACY': 'T', 'ACR': 'T', 'ACN': 'T', 
'AAT': 'N', 'AAC': 'N', 'AAY': 'N', 'AAA': 'K', 'AAG': 'K', 'AAR': 'K', 
'AGT': 'S', 'AGC': 'S', 'AGY': 'S', 'AGA': 'R', 'AGG': 'R', 'AGR': 'R', 
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GTY': 'V', 'GTR': 'V', 'GTN': 'V', 
'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GCY': 'A', 'GCR': 'A', 'GCN': 'A', 
'GAT': 'D', 'GAC': 'D', 'GAY': 'D', 'GAA': 'E', 'GAG': 'E', 'GAR': 'E', 
'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGY': 'G', 'GGR': 'G', 'GGN': 'G'
}

###################################################################################
def fasta_format(header, seq, lineLength = 60):
	'''format_fasta takes a header string and a sequence string. A '>' is
prepended to the header if not present. The sequence is broken into lines
of length lineLength (default 60). No testing is done for improper characters.
There should be no internal newlines in the header or the sequence.

This function returns a single string containing the header and all lines
of the sequence. The length of each line of the sequence is set to the
requested length.
'''
	if not header:
		header = '>Unnamed sequence'
	else:
		if not header.startswith('>'):
			header = '>' + header
		header = header.strip() #remove return if present
	
	temp = []
	temp.append(header)
	for x in xrange(0, len(seq), lineLength):
		temp.append(seq[x: x + lineLength])
	temp.append('')
	#add an extra, empty element so \n is added after last line
	fasta = '\n'.join(temp)
	return fasta
	#note: you can return as many references as you like
	#just separate them with commas

###################################################################################
def get_next_fasta (fileObject):
	'''usage: for header, seq in get_next_fasta(fileObject):
	
This is a generator that returns one fasta record's header and
sequence at a time from a multiple fasta file. Return character is removed
from the header. The sequence is returned as one continuous string
with no returns. The returned value is a tuple (header, sequence)
If their is no sequence associated with a header, seq will be an
empty string
Code simplification contributed by Dattatreya Mellacheruvu
01/16/2009, Jeffrey R. de Wet
08/02/2010 refactored to put lines into a temporary list
'''
	
	header = ''
	seq = ''

	lineList = []
	#The following for loop gets the header of the first fasta
	#record. Skips any leading junk in the file
	for line in fileObject:
		if line.startswith('>'):
			header = line.strip()
			break
	
	for line in fileObject:
		if line.startswith('>'):
			seq = ''.join(lineList)
			lineList = []
			yield header, seq
			header = line.strip()
			seq = ''
		else:
			#seq += line.strip()
			lineList.append(line.strip())
	#yield the last entry
	if header:
		seq = ''.join(lineList)
		yield header, seq


###################################################################################
def kmer_counter(seq, k):
	'''accepts a DNA sequence composed of ACGTNacgtn characters. Also accepts a parameter that is the length of the kmers to count. Always starts at
the beginning of a sequence and proceeds to the end, but does not include
any partial kmers. kmers overlap, in other words this is a sliding window
of length k that advances by on base at each step. Return a dictionary of
the kmers and the number of times each kmer was seen. The kmer keys in
the dictionary should be all uppercase, and when counting of the kmers
should not depend on the case of the sequence.
'''
	#define nsplit function that will split a string into a list of
	#smaller strings, each with k characters. the last element
	#has lefover characters
	
	def nsplit(seq, k):
		return [seq[m:m+k] for m in xrange(0, len(seq), k)]
	
	#use regular expression to check if sequence file is in a
	#acceptable form
	import re

	if not re.match(r'^[ACGTNacgtn]+$',seq):
		print("Failed! Your sequence does not match expectations!\n")
		return

	#define the dictionary of kmers and a list, used internally in function
	kmer_dictionary = dict()
	seqlist=[]
	
	#to have a sliding window, I create a list of sequences that will have
	#all possible reading frames to cover all possible kmers. The number of
	#sequences is dependent on the length of the kmer
	
	x=range(0,k)
	for i in x:
		seqlist.append(seq[i:])

	#Then I loop over all the reading frames, splitting the sequence into 		#a list of kmers (kmerslist)
	
	for seq in seqlist:

		kmerslist = nsplit(seq, k)
		
	#if the kmer is shorter that k (so called partial kmer), then 
	#it will not be tested for existance in dictionary.
	#otherwise all kmers will be added to dictionary and counted
	#testing for existance is case insensitive as all letters are
	#transformed into uppercase.
		
		for kmer in kmerslist:
			if len(kmer) < k:
				continue
			else:
				if kmer_dictionary.has_key(str.upper(kmer)):
					kmer_dictionary[str.upper(kmer)] += 1
				else:
					kmer_dictionary[str.upper(kmer)]=1
		
	return kmer_dictionary

################################################################################
def author_name ():
	'''returns my name, has no parameters
'''
	name= 'Kart Tomberg'
	return name


###################################################################################
def reverse_compliment(seq):
	'''accepts a string of ACGTNacgtn. The string may have a mix of upper and lower case characters, and the cases must be preserved in the returned sequence string.
'''

	#use regular expression to check if sequence file is in a
	#acceptable form
	import re

	if not re.match(r'^[ACGTNacgtn]+$',seq):
		print("Failed! Your sequence does not match expectations!\n")
		return

	#reverse the string (seq) given	
	reverse_seq = seq[::-1]

	#made a dictionary for complement nucleotides. I am using w,y instead 
	#of t,g as temporary replacements of the code to make sure that 
	#already replaced a's and c's would not be replaced back!
	complement = {
	'a':'w',
	'A':'W',
	'c':'y',
	'C':'Y',
	't':'a',
	'T':'A',
	'g':'c',
	'G':'C',
	'w':'t',
	'W':'T',
	'y':'g',
	'Y':'G',	
	}
	# I loop for sorted dictionary, to make sure a comes before t and c 
	#before g in order for my temporary replacement to work.

	for key in sorted(complement):
		reverse_seq=reverse_seq.replace(key,complement[key])

	return reverse_seq


###################################################################################
def base_counts(seq):
	'''accepts a seqeunce of ACGTNacgtn and returns the counts of As, Cs Gs Ts and Ns, regardless of there case in the input sequence. Return these as a dictionary with the counts indexed by the base.

'''

	#use regular expression to check if sequence file is in a acceptable 
	#form
	import re

	if not re.match(r'^[ACGTNacgtn]+$',seq):
		print("Failed! Your sequence does not match expectations!\n")
		return

	#turn string seq into a list where each nucleotide is an 
	#element	
	seq=list(seq)
	
	#define a dictionary called nucleotides
	nucleotides=dict()

	#counts different nucleotides (case insensitively) in the list	
	for nucleotide in seq:
		if nucleotides.has_key(str.upper(nucleotide)):
			nucleotides[str.upper(nucleotide)] += 1
		else:
			nucleotides[str.upper(nucleotide)]=1
		
	return nucleotides


###################################################################################
def gc_content(seq):
	'''accepts a DNA sequence composed of ACGTNacgtn and returns the percent GC content as a float.
The calculation of the percent does not take positions of variation (n,N) into account in the calculation.
float is rounded up with percision of 2 decimals and ranges from 0-100%.
'''
	#use regular expression to check if sequence file is in a 
	#acceptable form
	import re

	if not re.match(r'^[ACGTNacgtn]+$',seq):
		print("Failed! Your sequence does not match expectations!\n")
		return

	#use a dictionary called status_def to replace nucleotides in string 
	#with 0 or 1,based if are gc (=1) or not (=0). N status is defined as 
	#2 and is not taken into account in the calculations.
	
	status_def = {
	'a':'0',
	'A':'0',
	'c':'1',
	'C':'1',
	't':'0',
	'T':'0',
	'g':'1',
	'G':'1',
	'n':'2',
	'N':'2'	
	}
	
	#I make the replacement and return the sequence as a list called status
	for key in status_def:
		seq=seq.replace(key,status_def[key])
	status=list(seq)
	
	#count the number of gc-s and not gc-s
	no=status.count('0')
	yes=status.count('1')

	#calculate the gc percentage with the percision of 2 decimals	
	gc=round(float(yes)/float(yes+no)*100,2)
		
	return gc


###################################################################################
def translate_dna(seq,start_position=1,end_position=1):
	'''Accepts a sequence consisting of ACGTNacgtn and returns the encoded peptides in single letter code. Begins translating at the specified base position and ends at the specified base position (the first base will be considered position 1) Default is to start at the beginning of the sequence and translate to the end.
X will be used at any position where a codon is not in the table in addition to coding for stop.
'''
	#use regular expression to check if sequence file is in a
	#acceptable form
	import re

	if not re.match(r'^[ACGTNacgtn]+$',seq):
		print("Failed! Your sequence does not match expectations!\n")
		return

	#define x as the position in python because lists start with 0
	#if end is kept at default, y will be defined as length of the 
	#sequence, else the position given	
	x=start_position-1
	if end_position == 1:	
		y=len(seq)
	else:
		y=end_position

	#redefine seq based on the start and end position of translating
	seq=seq[x:y]

	#define nsplit function that will split a string into a list of smaller
	#strings, each with k characters. the last element has lefover 
	#characters
	
	def nsplit(seq, k):
		return [seq[m:m+k] for m in range(0, len(seq), k)]
	
	#use nsplit to make a list of codons (seqlist) of three nucleotides
	
	seqlist=nsplit(seq,3)

	#create a list called peptide and while looping over the seqlist, name 
	#the right aminoacid from the translation table and append it to the 
	#list called peptide. If not found in table, aminoacid is named 'X'
 
	peptide=[]
	for codon in seqlist:
		if transTab1L.__contains__(str.upper(codon)):
			aminoacid=transTab1L[str.upper(codon)]
		else:
			aminoacid='X'
		peptide.append(aminoacid)
	#join peptide into one string of aminoacids

	return "".join(peptide)	


###################################################################################
def reverse_translate(peptide):
	'''accepts an amino acid sequence and returns a dna sequence that encodes it. Because of the degeneracy of the code, you will need to use the IUPAC dna code to represent the degenerate positions. Returns dna sequence in lower case as a string.
'''
	#use regular expression to check if sequence file is in a acceptable 
	#form
	import re

	if not re.match(r'^[ABCDEFGHIKLMNPQRSTVWXYZ]+$',peptide):
		print("Failed! Your peptide does not match expectations!\n")
		return
	#turn the string into list of individual characters
	peptide=list(peptide)
	
	#made the reverse table for AA and its degenerative codons as a 
	#dictionary
	reverse_table = {
	'A': 'GCN',
	'B': 'RAY',
	'C': 'TGY',
	'D': 'GAY',
	'E': 'GAR',
	'F': 'TTY',
	'G': 'GGN',
	'H': 'CAY',
	'I': 'ATH',
	'K': 'AAR',
	'L': 'YTN',
	'M': 'ATG',
	'N': 'AAY',
	'P': 'CCN',
	'Q': 'CAR',
	'R': 'MGN',
	'S': 'TSN',
	'T': 'ACN',
	'V': 'GTN',
	'W': 'TGG',
	'X': 'NNN',
	'Y': 'TAY',
	'Z': 'SAR',
	}

	#define a list called seq and replace AA with its codons and append 
	#the codons into seq
	seq=[]
	for aminoacid in peptide:
		if reverse_table.has_key(aminoacid):
			codon=reverse_table[aminoacid]
			seq.append(codon)

	#join seq into a lowercase string
	return str.lower("".join(seq))

#!/usr/bin/python
import sys
import os
import fileinput
import itertools
from decimal import *
from datetime import datetime
import time
from math import sqrt

"""
TEAR allows the offline analysis of a sequence for repeat-induced point mutation using a profile based approach.
The user provides:
The sequences for use in the profile in fasta format (-p)
The sequence to be tested in fasta format (-i)
The file path for the output file (-o)
The 
"""

print 'Group Chi analysis and SD analysis combined'

# Needed variables

#profilefilePath = sys.argv[1]
profilefilePath = "/Users/Admin/Dropbox/Programs/tools/TEAR/temp/profile.fa"
#testInputSequence = sys.argv[2]
testInputSequence = "GCGATCGCGCGGCATCAGCATTATATTTTAAACCCGGGTTTAGCGATCATTATATACGGGCGCTCGTGTGTACACAC"
#outputPath = sys.argv[3]
outputPath = "/Users/Admin/Dropbox/Programs/tools/TEAR/temp/"
#theThreshold = sys.argv[4]
theThreshold = "1"

# setup any things needed later
bases = ['A', 'C', 'G', 'T']

def reverse(s):
	"""return the sequence string in reverse order"""
	letters = list(s)
	letters.reverse()
	return ''.join(letters)


def compliment(s):
	"""Return the complimentary sequence string"""
	basecompliment = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
	letters = list(s)
	letters = [basecompliment[base] for base in letters]
	return ''.join(letters)


def reversecompliment(s):
	"""return the reverse compliment of a DNA string."""
	s = reverse(s)
	s = compliment(s)
	return s

def gc(s):
	"""return the percentage of the DNA string composed of GC."""
	gc = s.count('G') + s.count('C')
	return gc * 100.0 / len(s)
	
def kmer_count(s):
	"""Creates a dictionary containing the occurrence of k-mers of a sequence"""
	x = 0
	usrin = dict()
	while x != len(s):
		kmer  = s[(0+x):(int(kval)+x)]
		if kmer not in usrin:
			usrin[kmer] = 1
		else:
			usrin[kmer] += 1
		x = x+1
	for k,v in usrin.items():
		km = ((100.0/len(s))*v)
		usrin[k] = km
	for i, m in usrin.items():
		if len(i) == 1:
			del usrin[i]
	return usrin
	
def round_to_N(s, n):
	n = int(n)
	str_version = str(s)
	rounded = str_version[0:n]
	return float(rounded)

def input_parser(file_path):
	# open as fasta and remove illegal chars
	input_file = open(file_path, "r")
	
	fasta_dict = {}
	for line in input_file:
		if line[0:1] == ">":
			title = line
			title = title.replace("\n","")
			title = title.replace(">","")
			fasta_dict[title] = ""
		else:
			seq = line
			seq = seq.strip()
			fasta_dict[title] = fasta_dict[title] + seq
	return fasta_dict
			



# setting the scene	

kval = 2 # The length of k-mer
SD_value = float(theThreshold)


#  1- getting the sequences to be used for the profile

tempfile = input_parser(profilefilePath)

tempnewlist = tempfile
print tempnewlist

table = {}

j = 1

for key in tempnewlist:
	seq_name = 'seq_%s' % (j)
	pre_seq = tempnewlist[key]
	post_seq = pre_seq.upper()
	for char in ['N','W','K','M','S','X','R','Y','D','H','P','Q','U','I','O','V','B','Z']:
		if char in post_seq:
			post_seq = post_seq.replace(char, current_placeholder)
	table[seq_name] = post_seq
	j = j + 1

print "The input sequences are: "
print ""
print table
print ""
print "length of table "
print len(table)

# The test sequence and remove illegal chars

test_input_pre = str(testInputSequence)
tempthis_N = test_input_pre.upper()
for char in ['N','W','K','M','S','X','R','Y','D','H','P','Q','U','I','O','V','B','Z']:
	if char in tempthis_N:
		tempthis_N = tempthis_N.replace(char, current_placeholder)
test_input = tempthis_N


#  2- count different kmers and add results to matrix

AA = []
AT = []
AG = []
AC = []
TT = []
TA = []
TG = []
TC = []
GA = []
GT = []
GC = []
GG = []
CA = []
CT = []
CC = []
CG = []
#for bugs
A = []
C = []
T = []
G = []

# this table of lists is where the results are stored
kmer_boxdict = {'AA':[], 'AT':[], 'AG':[], 'AC':[], 'TA':[], 'TT':[], 'TG':[], 'TC':[], 'GA':[], 'GT':[], 'GC':[], 'GG':[], 'CA':[], 'CT':[], 'CG':[], 'CC':[], 'A':[], 'T':[], 'G':[], 'C':[]}



# this takes each item in the table of inputs and places it and its respective value in the kmer_boxdict 
print "---------------------The kmer counts for the table --------------------------"

for k, v in table.items():
	print k, v
	print ""
	test1 = kmer_count(table[k])
	for r, t in test1.items():
		#print r
		kmer_boxdict[r].append(t)

print "---------------------The kmer counts for the input seq --------------------------"
input_dict = kmer_count(test_input)



#test1 = kmer_count(table['seq_2'])
#print test1

for k, v in input_dict.items():
	print k, v

print ""
print ""
print "kmer_boxdict "
for i, o in kmer_boxdict.items():
	print i, o
	print ""

# 3- conduct statistics
# calculate SD for each kmer




print "-------------------------------------- The variance is: ----------------------------------"

# this table of lists is where the results are stored for the variance
kmer_var_dict_values = {'AA':[], 'AT':[], 'AG':[], 'AC':[], 'TA':[], 'TT':[], 'TG':[], 'TC':[], 'GA':[], 'GT':[], 'GC':[], 'GG':[], 'CA':[], 'CT':[], 'CG':[], 'CC':[], 'A':[], 'T':[], 'G':[], 'C':[]}

# this table of lists is where the results are stored for the variance
kmer_var_dict = {'AA':[], 'AT':[], 'AG':[], 'AC':[], 'TA':[], 'TT':[], 'TG':[], 'TC':[], 'GA':[], 'GT':[], 'GC':[], 'GG':[], 'CA':[], 'CT':[], 'CG':[], 'CC':[], 'A':[], 'T':[], 'G':[], 'C':[]}


for i, o in kmer_boxdict.items():
	print i, o
	denom = len(o)
	sumel = sum(o)
	print denom
	print sumel
	if denom != 0:
		sam_mean = sumel/denom
		print " mean: "
		print round_to_N(sam_mean,8)
		for j in o:
			print j
			diff = (round_to_N(j,6) - round_to_N(sam_mean,6))**2
			print diff
			kmer_var_dict_values[i].append(diff)

print "---------------------- the dict of the variance values (squared diff from the mean) looks like so ------------------------------------"
 

for i , o in kmer_var_dict_values.items():
	print i
	print o

for i , o in kmer_boxdict.items():
	print "i = "
	print i
	print "o = "
	print o
	if len(o) != 0:
		kmer_var_dict[i].append((sum(kmer_var_dict_values[i])/(len(kmer_var_dict_values[i])-1)))
	else:
		print "=0"


print "---------------------- the dict of the variance looks like so ------------------------------------"
for i, o in kmer_var_dict.items():
	print i, o
# compute the average of the variables and take the square root


# this table of lists is where the results are stored
kmer_SD_dict = {'AA':[], 'AT':[], 'AG':[], 'AC':[], 'TA':[], 'TT':[], 'TG':[], 'TC':[], 'GA':[], 'GT':[], 'GC':[], 'GG':[], 'CA':[], 'CT':[], 'CG':[], 'CC':[], 'A':[], 'T':[], 'G':[], 'C':[]}


for i, o in kmer_var_dict.items():
	if len(o) != 0:
		tempvalue = kmer_var_dict[i]
		print tempvalue[0]
		diff2 = sqrt(tempvalue[0])
		kmer_SD_dict[i].append(diff2)
		
print "--------------------------------------- The SD is: ------------------------------ "
for i, o in kmer_SD_dict.items():
	print i, o

# compare query to profile
#print "The input sequence"
#for i, o in input_dict.items():
	#print i, o

# dictionaries:
#Profile sequence composition:			kmer_var_dict
#For mean do:							sum(kmer_var_dict)/len(kmer_var_dict)
#SD of the profile sequences:			kmer_SD_dict
#input sequence composition:			input_dict


results_dict = {'AA':[], 'AT':[], 'AG':[], 'AC':[], 'TA':[], 'TT':[], 'TG':[], 'TC':[], 'GA':[], 'GT':[], 'GC':[], 'GG':[], 'CA':[], 'CT':[], 'CG':[], 'CC':[], 'A':[], 'T':[], 'G':[], 'C':[]}

#Variables:
#The threshold
the_SD = SD_value
#print "the test"
#print "With %s SD applied" % (the_SD)
for g, h in input_dict.items():
	if len(kmer_SD_dict[g]) != 0:
		local_SD = kmer_SD_dict[g]
		local_mean = sum(kmer_boxdict[g])/len(table)
		if ( ( h > (local_mean + local_SD[0]*the_SD)) or (h < (local_mean - local_SD[0]*the_SD))  ):
			print "Query value = "
			print g, h
			print "Profile value = "
			print kmer_var_dict[g]
			print (sum(kmer_var_dict[g]))/len(kmer_var_dict[g])
			print "SD"
			print kmer_SD_dict[g]
			results_dict[g].append("Fail")
			print " "
		else:
			print "%s value is within SD" % g
			results_dict[g].append("Pass")
			print " "
			
			
# report kmers with significant differences
print "results:"

phpResultsList = ""

#for i, o in input_dict.items():
#	print i, o

#for i, o in results_dict.items():
#	print i, o

print "-------------------------Chi squared analysis--------------------------- "

Chi_dict = {'AA':[], 'AT':[], 'AG':[], 'AC':[], 'TA':[], 'TT':[], 'TG':[], 'TC':[], 'GA':[], 'GT':[], 'GC':[], 'GG':[], 'CA':[], 'CT':[], 'CG':[], 'CC':[], 'A':[], 'T':[], 'G':[], 'C':[]}

for i, o in input_dict.items():
	obs = input_dict[i]
	obs = round_to_N(obs,6)
	exp = (sum(kmer_boxdict[i])/len(table))
	exp = round_to_N(exp,6)
	chi_value = (((obs - exp)**2)/4)
	Chi_dict[i].append(chi_value)

print "-------------------------Analysis complete, formatting--------------------------- "

print ""

print ""


for i, o in input_dict.items():
	if len(kmer_boxdict[i]) != 0:
		Mean = (sum(kmer_boxdict[i])/len(table))
	else:
		Mean = "0"
	SD = kmer_SD_dict[i]
	isSignificant = results_dict[i]
	Chi = Chi_dict[i]
	print phpResultsList
	phpResultsList = phpResultsList + "%s,%s,%s,%s,%s," % (i, o, Mean, SD[0], isSignificant[0])
print phpResultsList
print "final"

unique_jod_id = str(unique_jod_id)

save_file_path = outputPath + unique_jod_id + '.txt'
temp_out_file = open(save_file_path, 'w')
temp_out_file.write(phpResultsList)

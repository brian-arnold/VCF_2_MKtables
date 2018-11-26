#!/usr/bin/python -tt

""" 
INSERT USEFUL COMMENTS HERE
"""

import sys
import vcf 
import re

def getChromList(bedList):
	ChromList = []
	# Only go through first file, as I assume all BED files will have all chromosomes represented
	bedFH = open(bedList[0], 'r')
	for line in bedFH:
		l = line.split()
		chrom = l[0]
		if chrom not in ChromList:
			ChromList.append(chrom)
	return(ChromList)

def getBedIntervalsDict(bedList):
	# Open bed files
	# NOTE: bedfiles have 0-based coordinates, half open, meaning it doesn't include ending coordinate
	# convert to 1-based closed interval by adding 1 to start position and taking rest of coorinates as is
	bedIntervalsDict = {}
	for bedFN in bedList:
		bedIntervalsDict[bedFN] = {}
		#print(bedFN)
		bedFH = open(bedFN, 'r')
		for line in bedFH:
			l = line.split()
			chrom = l[0]
			start = int( l[1] ) + 1 
			end = int( l[2] )
			try:
				bedIntervalsDict[bedFN][chrom].update(range(start, end+1))
			except:
				bedIntervalsDict[bedFN][chrom] = set(range(start, end+1))	# range doesn't include endpoint, so +1

			#print(chrom, start, end)
		bedFH.close()
	return(bedIntervalsDict)

def getCallableSites(chrList, bedList, bedIntDict):
	CallableSitesDict = {}
	for chrom in chrList:
		CallableSitesDict[chrom] = {}
		first = 1
		TempSet = set()
		for bedFN in bedList:
			if first==1:
				TempSet = set(bedIntDict[bedFN][chrom])
				first = 0
				#print(bedIntDict[bedFN][chrom])
				#print(bedFN, chrom, TempSet)
			else:
				TempSet = TempSet.intersection(set(bedIntDict[bedFN][chrom]))
				#print(bedIntDict[bedFN][chrom])
				#print(bedFN, chrom, TempSet)
		
		#print(TempSet)
		for x in TempSet:
			CallableSitesDict[chrom][x]=1	

	return(CallableSitesDict)

def main():
	# initialization

	# parse command-line arguments
	if len(sys.argv) == 3:
		# have comma separated lists for each field, split on comma
		vcfTempList = sys.argv[1]
		bedTempList = sys.argv[2]

		vcfList = vcfTempList.split(",")
		bedList = bedTempList.split(",")
	else:
		print("USAGE: Specify two things, comma-separated (no spaces) list of vcfs, comma-separated (no spaces) list of bed files")
		sys.exit()
	
	# get list of chromosomes
	ChromList = getChromList(bedList) # list of chromosome names

	# Open bed files
	# NOTE: bedfiles have 0-based coordinates, half open, meaning it doesn't include ending coordinate
	# convert to 1-based closed interval by adding 1 to start position and taking rest of coorinates as is
	bedIntervalsDict = getBedIntervalsDict(bedList) # 2D Dict of sets, bedIntervalsDict[bedFN][chrom] = set(range(start, end+1))

	# for each chrom, find intersection of all intervals across all input bed files
	# load into 2D dict CallableSitesDict, CallableSitesDict[chrom][pos]=1
	CallableSitesDict = getCallableSites(ChromList, bedList, bedIntervalsDict) 

	RefDict = {}	# 2D dict, RefDict[chromosome][site] = reference allele
	AltDict = {}	# 3D dict, AltDict[vcfFile][chromosome][site] = (alternate allele, frequency)
	MultiDict = {} 	# Dict recording multiallelic sites and indels
	# Initialize Dicts
	for chrom in CallableSitesDict:
		MultiDict[chrom] = {}
		RefDict[chrom] = {}
		for vcfFN in vcfList:
			AltDict[vcfFN] = {}
			AltDict[vcfFN][chrom] = {}

	# Go through each VCF file
	for vcfFN in vcfList:	
		myVcf = vcf.Reader(filename=vcfFN)
		# Go through ea line in VCF
		for record in myVcf:
			if record.CHROM in CallableSitesDict:
				if record.POS in CallableSitesDict[record.CHROM]:
					ref = str(record.REF)
					if re.match(r'[AGCT]$', ref): # ensure REF is a single base; indels  could make it longer string
						RefDict[record.CHROM][record.POS] = ref # occasionally gets overwritten
						if record.FILTER == []:	# it seems PyVCF sets FILTER to empty if it's PASS
							alt = str(record.ALT)
							if re.match(r'\[[AGCT]{1}\]', alt): # biallelic? ALT must match only 1 of A,G,C,T
								AltDict[vcfFN][record.CHROM][record.POS] = (alt, record.INFO['AF'])
							else:
								MultiDict[record.CHROM][record.POS] = 1	

	# Print output	
	outFile = open('InfoForMKtable.txt', 'w')
	print("Chromosome", "Position", "RefAllele", end=" ", file=outFile)
	for vcfFN in sorted(vcfList):
		print(vcfFN, "_Allele ", vcfFN, "_Freq",sep="", end=" ", file=outFile)
	print(file=outFile)
	for chrom in RefDict:
		for pos in RefDict[chrom]:
			# ignore multiallelic sites and indels
			if pos not in MultiDict[chrom]:
				print(chrom, pos, end=" ", file=outFile)
				print(RefDict[chrom][pos], end=" ", file=outFile)
				for vcfFN in sorted(vcfList):
					if pos in AltDict[vcfFN][chrom]:
						# remove square brackets PyVCF puts on ALT and AF
						x = re.sub('[\[\]]', '', str(AltDict[vcfFN][chrom][pos][0]))
						y = re.sub('[\[\]]', '', str(AltDict[vcfFN][chrom][pos][1]))
						print(x, y, file=outFile)
					else:
						print(RefDict[chrom][pos], 1.0, file=outFile)

if __name__ == '__main__':
  main()

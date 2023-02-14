#!/usr/bin/env python3

import sys
import re
import argparse
import gzip
import matplotlib.pyplot as plt

def argparser():
	parser = argparse.ArgumentParser()
	parser.add_argument("--methylkit", required = True, help = "gzipped .methylKit file")
	parser.add_argument("--fasta", required = False, help = "Reference fasta file")
	parser.add_argument("--binary_fasta", required = False, help = "'Binary' fasta file (CpGs as 1, rest as 0)")
	parser.add_argument("--binary_methylkit", required = False, help = "'Binary' fasta with methylation positions")
	parser.add_argument("--all_bins", required = False, help = "File with ")
	parser.add_argument("--prefix", required = True, help = "Prefix for the output file")
	parser.add_argument("--annotation", required = True, help = "Annotation GTF file")
	parser.add_argument("--plot", required = False, help = "Whether you want to plot the results or not")
	
	args = parser.parse_args()
	return args

def write_fasta(fasta_dict, file_handle, string):
	''' Writes a dict of sequences into a fasta file'''
	
	print("Writing {}".format(file_handle))
	out = open(file_handle, "w")
	for chrom in fasta_dict:
		if string == True:
			out.write(">{}\n{}\n".format(chrom, fasta_dict[chrom]))
		else:
			out.write(">{}\n{}\n".format(chrom, "".join(map(str, fasta_dict[chrom]))))
			

def parse_fasta(fasta):
	''' Parses a fasta file and stores it in a dictionary'''

	print("Parsing the fasta file...")
	fasta_dict = {}
	if fasta.endswith(".gz"):
		fasta = gzip.open(fasta, "rb")
		for line in fasta:
			line = str(line, "utf-8")
			if line.startswith(">"):
				chrom = line.strip().split(" ")[0][1:]
				print(chrom)
				fasta_dict[chrom] = []
			else:
				fasta_dict[chrom].append(line.strip())
		fasta.close()
	else:
		fasta = open(fasta, "r")
		for line in fasta:
			if line.startswith(">"):
				chrom = line.strip().split(" ")[0][1:]
				print(chrom)
				fasta_dict[chrom] = []
			else:
				fasta_dict[chrom].append(line.strip())
		fasta.close()
	return fasta_dict

def convert_to_binary(fasta_dict):
	''' 
	Reads a fasta dictionary and converts it to 'binary', in which
	every CG position is a 1, and all other bases are 0s
	'''

	print("Converting the fasta to 'binary'")
	bin_dict = {}
	cpg_bin_dict = {}
	for chrom in fasta_dict:
		print(chrom)
		seq = "".join(fasta_dict[chrom])
		bin_dict[chrom] = [0] * len(seq)
		cpg_bin_dict[chrom] = [0] * len(seq)
		cgs = [cg.start() for cg in re.finditer("CG", seq)]
		for pos in cgs:
			bin_dict[chrom][pos] = 1
			bin_dict[chrom][pos + 1] = 1
	
	return bin_dict, cpg_bin_dict
		
def parse_binary_fasta(binary_fasta):
	''' Parses a pre-made 'binary' fasta, assuming the sequences are on one line each'''

	print("Parsing the 'binary' fasta file...")
	bin_dict = {}
	cpg_bin_dict = {}
	if binary_fasta.endswith(".gz"):
		binary_fasta = gzip.open(binary_fasta, "rb")
		for line in binary_fasta:
			line = str(line, "utf-8")
			if line.startswith(">"):
				chrom = line.strip().split(" ")[0][1:]
				bin_dict[chrom] = ""
			else:
				bin_dict[chrom] = line.strip()
		binary_fasta.close()

	else:
		binary_fasta = open(binary_fasta, "r")
		for line in binary_fasta:
			if line.startswith(">"):
				chrom = line.strip().split(" ")[0][1:]
				bin_dict[chrom] = ""
			else:
				bin_dict[chrom] = line.strip()
		binary_fasta.close()
	
	return bin_dict, cpg_bin_dict

def make_cpg_dict(bin_dict):
	''' 
	'''

	print("Making the CpG dict")
	cpg_bin_dict = {}
	for chrom in bin_dict:
		print(chrom)
		cpg_bin_dict[chrom] = [0] * len(bin_dict[chrom])

	return cpg_bin_dict

def parse_all_bins(all_bins_file):
	
	print("Parsing the bins file")
	file = open(all_bins_file, "r")
	for line in file:
		all_bins = line.strip().split("\t")
	
	all_bins = [int(x) for x in all_bins]
	
	return(all_bins)
	


def parse_annotation(annotation):
	''' Parses annotation gtf file to find coordinates of TSSs +/- 2kb'''

	print("Parsing the annotation file...")
	tss_dict = {}
	annotation = open(annotation, "r")
	for line in annotation:
		if not line.startswith("##"):
			line = line.split("\t")
			if line[2] == "transcript":
				chrom = line[0]
				if chrom not in tss_dict:
					print(chrom)
					tss_dict[chrom] = []
				if line[6] == "+":
					start = int(line[3]) - 1
					start -= 2000
					end = start + 4000
					tss_dict[chrom].append(((start, end), line[6]))
				else:
					start = int(line[4])
					start -= 2000
					end = start + 4000
					tss_dict[chrom].append(((start, end), line[6]))
	annotation.close()
	
	return tss_dict	

def parse_methylkit(methylkit, cpg_bin_dict):
	''' Turns any 0s in cpg_bin_dict to 1s if present in methylkit file'''

	print("Parsing the methylkit file...")
	if methylkit.endswith(".gz"):
		methylkit_file = gzip.open(methylkit, "rb")
	else:
		methylkit_file = open(methylkit, "r")
	for line in methylkit_file:
		if methylkit.endswith(".gz"):
			line = str(line, "utf-8")		
		if not line.startswith("chrBase"):
			line  = line.split("\t")
			pos = int(line[2]) - 1
			chrom = line[1]
			cpg_bin_dict[chrom][pos] = 1
	
	for chrom in cpg_bin_dict:
		temp = "".join(map(str, cpg_bin_dict[chrom]))
		cpg_bin_dict[chrom] = temp
		
	return cpg_bin_dict
		

def calculate_meth_bin(tss_dict, cpg_bin_dict):
	'''
	Checks the bins around TSSs, if there is coverage in the bin then
	adds that bin to the methylation_bins dictionary, keeping track of the
	total number of bins with coverage
	'''

	methylation_bins = [0 for x in range(400)]
	for chrom in tss_dict:
		print(chrom)
		for pos in tss_dict[chrom]:
			meth_slice = cpg_bin_dict[chrom][pos[0][0]:pos[0][1]]
			my_bin = 0
			for x in range(0,len(meth_slice), 10):
# 				if sum(meth_slice[x:x+10]) != 0:
				if "1" in meth_slice[x:x+10]:
					if pos[1] == "+":
						methylation_bins[my_bin] += 1
					else:
						methylation_bins[399 - my_bin] += 1
				my_bin += 1
	
	return methylation_bins
	
def calculate_both_bins(tss_dict, bin_dict, cpg_bin_dict):
	''' Will calculate both all_bins and methylation_bins if neither is specified'''

	print("Calculating the bins...")
	all_bins = [0] * 400
	methylation_bins = [0] * 400
	for chrom in tss_dict:
		print(chrom)
		for pos in tss_dict[chrom]:
			my_slice = bin_dict[chrom][pos[0][0]:pos[0][1]]
			meth_slice = cpg_bin_dict[chrom][pos[0][0]:pos[0][1]]
			my_bin = 0
			for x in range(0,len(my_slice), 10):
				if "1" in my_slice[x:x+10]:
					if pos[1] == "+":
						all_bins[my_bin] += 1
					else:
						all_bins[399 - my_bin] += 1
				my_bin += 1
		
			my_bin = 0
			for x in range(0,len(meth_slice), 10):
				if "1" in meth_slice[x:x+10]:
					if pos[1] == "+":
						methylation_bins[my_bin] += 1
					else:
						methylation_bins[399 - my_bin] += 1
				my_bin += 1
	
	return all_bins, methylation_bins

def normalize_bins(all_bins, methylation_bins):
	''' Normalizes the methylation bin to percent of bins covered that have CpGs'''
	
	print("Normalizing the bins")	
	normal_all_bins = [(all_bins[x] / all_bins[x]) * 100 for x in range(len(all_bins))]
	normal_methylation_bins = [(methylation_bins[x] / all_bins[x]) * 100 for x in range(len(all_bins))]

	return(normal_all_bins, normal_methylation_bins)	

def plot_bins(all_bins, methylation_bins):
	''' Plots the bins simply'''
	
	print("Plotting!")
	plt.plot(all_bins)
	plt.plot(methylation_bins)
	plt.show()

def write_files(all_bins, methylation_bins, normal_methylation_bins, prefix):
	''' Writes the bins to tab delimited files'''

	print("Writing files...")

	out = open("{}_allbins_raw.tab".format(prefix), "w")
	out.write("\t".join(map(str, all_bins)))
	out.close()
	
	out = open("{}_methylation_raw.tab".format(prefix), "w")
	out.write("\t".join(map(str, methylation_bins)))
	out.close()
	
	out = open("{}_methylation_normalized.tab".format(prefix), "w")
	out.write("\t".join(map(str, normal_methylation_bins)))
	out.close()
	
	

def run_script():
	print("Running script...")
	args = argparser()
	
	if args.fasta == None and args.binary_fasta == None:
		print("Need either a fasta or the 'binary' fasta file!")
		sys.exit()
	
	if args.binary_fasta:
		bin_dict = parse_fasta(args.binary_fasta)
		if args.binary_methylkit:
			cpg_bin_dict = parse_fasta(args.binary_methylkit)
		else:
			cpg_bin_dict = make_cpg_dict(bin_dict)
			cpg_bin_dict = parse_methylkit(args.methylkit, cpg_bin_dict)
	elif args.fasta:
		fasta_dict = parse_fasta(args.fasta)
		bin_dict, cpg_bin_dict = convert_to_binary(fasta_dict)
		write_fasta(bin_dict, "{}_cpgs.fa".format(args.prefix), False)
	
	tss_dict = parse_annotation(args.annotation)
	write_fasta(cpg_bin_dict, "{}_covered_cpgs.fa".format(args.prefix), False)
	
	if args.all_bins:
		all_bins = parse_all_bins(args.all_bins)
		methylation_bins = calculate_meth_bin(tss_dict, cpg_bin_dict)
	else:
		all_bins, methylation_bins = calculate_both_bins(tss_dict, bin_dict, cpg_bin_dict)
	
	normal_all_bins, normal_methylation_bins = normalize_bins(all_bins, methylation_bins)
	
	write_files(all_bins, methylation_bins, normal_methylation_bins, args.prefix)
	
	if args.plot:
		plot_bins(all_bins, methylation_bins)
		plot_bins(normal_all_bins, normal_methylation_bins)
	
if __name__ == "__main__":
	run_script()
	print("All done!")
	



#!/usr/bin/python

########################################################################
# 11 Oct 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################


import sys, re, os
import subprocess
import argparse
import pkg_resources

def read_tfs(ifile):
	data = {}
	with open(ifile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			data[word[0]] = 1
	return data

def select_tfs(ifile, ofile, tfs):
	output = open(ofile, "w")
	with open(ifile) as f:
		header= next(f)
		header = header.rstrip()
		output.write("{}\n".format(header)),
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[0] in tfs:
				output.write("{}\n".format(line)),

def main():
	parser = argparse.ArgumentParser(description='Select mouse TF\'s from a expression table\n')
	parser.add_argument('-i','--input', help='Input file in tab delimiated format containing ensembl IDs on the first column', required=True)
	parser.add_argument('-o','--outfile', help='Output file', required=True)
	args = vars(parser.parse_args())

	tfs = pkg_resources.resource_filename('pygprofile', 'data/animal_tfdb_mouse.txt')
	tfs_dict = read_tfs(tfs)
	select_tfs(args["input"], args["outfile"], tfs_dict)
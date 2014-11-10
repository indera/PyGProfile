#!/usr/bin/python

########################################################################
# 28 July 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################


import sys, re, os
import subprocess
import argparse
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
#import pygenerich
import rpy2.robjects.numpy2ri as rpyn
import numpy as np
import pkg_resources
import pygenerich
from pygenerich.wrappers import spia, custom_pathway_analysis

def _dict_to_namedvector(init_dict):
	"""Call R to create a named vector from an input dictionary."""
	return ro.r.c(**init_dict)

def read_input_file(input_file, lfc=False):
	gene_vals = {}
	with open(input_file) as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[1]=="NA":
				pass
			else:
				if lfc:
					gene_vals[word[0]] = (word[1], word[2])
				else:
					gene_vals[word[0]] = float(word[2])
	return gene_vals

def read_deseq_file(input_file, lfc=False, adjpval=False):
	gene_vals = {}
	with open(input_file) as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[5]=="NA":
				pass
			else:
				if adjpval:
					if lfc:
						gene_vals[word[0]] = (word[2], word[6])
					else:
						gene_vals[word[0]] = float(word[6])
				else:
					if lfc:
						gene_vals[word[0]] = (word[2], word[6])
					else:
						gene_vals[word[0]] = float(word[6])
	return gene_vals

def read_gfold_file(input_file, lfc=False):
	#Unsure what to do here since the result depends on gfold and not p-value. 
	gene_vals = {}
	with open(input_file) as f:
		next(f)
		for line in f:
			if line.startswith("#"):
				pass
			line = line.rstrip()
			word = line.split("\t")
			if lfc:
				gene_vals[word[0]] = (0, word[3])
			else:
				gene_vals[word[0]] = float(word[2])
	return gene_vals

def main():
	parser = argparse.ArgumentParser(description='Pathway enrichment analysis script. Please choose which program you want to run. Custom input file must bed in the format "IDs\tLFC\tP-value"\n')
	parser.add_argument('-i','--input', help='File containing ID\tLFC\tP-value', required=False)
	parser.add_argument('-d','--deseq', help='DESEQ2 output file', required=False)
	#parser.add_argument('-g','--gfold', help='Gfold Input file. Not recommended!', required=False)
	parser.add_argument('-p','--pval', help='P-value filter for significance', required=False)
	parser.add_argument('-a','--apval', help='Adjusted P-value filter for significance, DESEQ2 only', required=False)
	parser.add_argument('-spia',action='store_true', help='SPIA, only use if you are comparing conditions!',required=False)
	parser.add_argument('-path', action='store_true', help='Hypergeometric test for pathways', required=False)
	args = vars(parser.parse_args())
	
	path = os.path.dirname(gene_enrichment_toolkit.__file__)
	kegg_info = pkg_resources.resource_filename('pygprofile', 'data/mouse_kegg.gmt')
	if args["spia"]:
		if args["input"]:
			gene_vals = read_input_file(args["input"], lfc=True)
			spia.spia(gene_vals, args["pval"])
		elif args["deseq"]:
			if args["pval"]:
				gene_vals = read_deseq_file(args["deseq"], lfc=True)
				spia.spia(gene_vals, args["pval"])
			elif args["apval"]:
				gene_vals = read_deseq_file(args["deseq"], lfc=True, adjpval=True)
				spia.spia(gene_vals, args["apval"])
	if args["path"]:
		if args["input"]:
			gene_vals = read_input_file(args["input"], lfc=False)
			custom_pathways_analysis.custom_pathways(gene_vals,kegg_info, args["pval"])
		elif args["deseq"]:
			if args["pval"]:
				gene_vals = read_deseq_file(args["deseq"])
				custom_pathways_analysis.custom_pathways(gene_vals,kegg_info, args["pval"])
			elif args["apval"]:
				gene_vals = read_deseq_file(args["deseq"], adjpval=True)
				custom_pathways_analysis.custom_pathways(gene_vals,kegg_info, args["apval"])

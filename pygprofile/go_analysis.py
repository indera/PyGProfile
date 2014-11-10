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
import rpy2.robjects.numpy2ri as rpyn
import numpy as np
from collections import defaultdict
from pygenerich.wrappers import topgo, gostats, goprofiles

def _dict_to_namedvector(init_dict):
	"""Call R to create a named vector from an input dictionary."""
	return ro.r.c(**init_dict)

def read_input_file(input_file):
	gene_vals = {}
	c = 0
	with open(input_file) as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 2:
				if word[1]=="NA":
					pass
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


def filter_genes_by_go(ifile, go_to_search, out):
	data = {}
	#BiomaRt working:
	go_ids = defaultdict(list)
	genome="mmusculus_gene_ensembl"
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	filter1 = "ensembl_gene_id"
	with open(ifile) as f:
#		header = next(f)
#		print header,
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			data[word[0]] = line
	values = []
	for key in data:
		values.append(key)
	#Need to make sure by looking at webservice if the numbers of genes are the same!
	#Would also be good to include name of GO term!
	#Also genes can have multiple GO terms so might use a list for dict!
	goooo = ro.r.getBM(attributes=StrVector([filter1, "go_id", "external_gene_name"]), filters="go_parent_term", values=go_to_search, mart=ensembl) #Need to use parent because go_id does not give everything!
	genes = goooo.rx(True,1)
	gos = goooo.rx(True,2)
	ext = goooo.rx(True,3)
	for x in range(0, len(genes)):
		go_ids[genes[x]] = (gos[x], ext[x])
	output = open(out, "w")
	for key in go_ids:
		output.write("{}\t{}\t{}\n".format(data[key], go_ids[key][0], go_ids[key][1])),
	output.close()

def main():
	parser = argparse.ArgumentParser(description='Gene Ontology and Pathway enrichment analysis script. Please choose which program you want to run. Custom input file must bed in the format "IDs\tLFC\tP-value"\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	go_parser = subparsers.add_parser('go', help="Do GO analysis on list of genes")
	search_parser = subparsers.add_parser('search', help="Search for GO genes")
	go_parser.add_argument('-i','--input', help='File containing ID\tP-value', required=False)
	go_parser.add_argument('-d','--deseq', help='DESEQ2 output file', required=False)
	go_parser.add_argument('-g','--gfold', help='Gfold Input file. Not recommended!', required=False)
	go_parser.add_argument('-p','--pval', help='P-value filter for significance', required=False)
	go_parser.add_argument('-a','--apval', help='Adjusted P-value filter for significance, DESEQ2 only', required=False)
	go_parser.add_argument('-topgo',action='store_true', help='TopGO', required=False)
	go_parser.add_argument('-gostats',action='store_true', help='GOstats', required=False)
	go_parser.add_argument('-goprofiles', action='store_true', help='goProfiles', required=False)
	go_parser.add_argument('-outname', help='Output name for results', required=False)
	search_parser.add_argument('-i','--input', help='Input matrix, ensembl ids must be on first column')
	search_parser.add_argument('-g', '--go', help='GO terms to search for, can be multiple', nargs="+")
	search_parser.add_argument('-o', '--output', help='Output file name', required=True)
	args = vars(parser.parse_args())
	if args["subparser_name"] == "go":
		if args["topgo"]:
			if args["input"]:
				gene_vals = read_input_file(args["input"])
				topgo.run_topgo(gene_vals, args["pval"], args["outname"])
			elif args["deseq"]:
				if args["apval"]:
					gene_vals = read_deseq_file(args["deseq"], adjpval=True)
					topgo.run_topgo(gene_vals, args["apval"], args["outname"])
				else:
					gene_vals = read_deseq_file(args["deseq"])
					topgo.run_topgo(gene_vals, args["pval"], args["outname"])
		if args["gostats"]:
			if args["input"]:
				gene_vals = read_input_file(args["input"])
				gostats.go_stats(gene_vals, args["pval"], args["outname"])
			elif args["deseq"]:
				if args["apval"]:
					gene_vals = read_deseq_file(args["deseq"], adjpval=True)
					gostats.go_stats(gene_vals, args["apval"], args["outname"])
				else:
					gene_vals = read_deseq_file(args["deseq"])
					gostats.go_stats(gene_vals, args["pval"], args["outname"])
		if args["goprofiles"]:
			if args["input"]:
				gene_vals = read_input_file(args["input"])
				goprofiles.goprofiles(gene_vals, args["pval"], args["outname"])
			elif args["deseq"]:
				if args["apval"]:
					gene_vals = read_deseq_file(args["deseq"], adjpval=True)
					goprofiles.goprofiles(gene_vals, args["apval"], args["outname"])
				else:
					gene_vals = read_deseq_file(args["deseq"])
					goprofiles.goprofiles(gene_vals, args["pval"], args["outname"])
	elif args["subparser_name"] == "search":
		filter_genes_by_go(args["input"], args["go"], args["output"])

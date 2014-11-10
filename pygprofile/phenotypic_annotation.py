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
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector

def read_input(ifile):
	genes = {}
	with open(ifile) as f:
		header = next(f)
		header = header.rstrip()
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			genes[word[0]] = line
	return genes, header

def annotate_ensembl(dict_obj):
	#Genes must be the keys
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	genome="mmusculus_gene_ensembl"
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	values = []
	for key1 in dict_obj:
		values.append(key1)
	C1BM = ro.r.getBM(attributes=StrVector(["ensembl_gene_id", "mgi_id"]), filters="ensembl_gene_id", values=values, mart=ensembl)
	gene = list(C1BM.rx(True,1))
	ent = list(C1BM.rx(True,2))
	data = {}
	for index, g in enumerate(gene):
		data[g] = ent[index]
	return data

def read_gene_anno(ifile):
	g_anno = defaultdict(list)
	with open(ifile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			mgi = word[5].strip()
			mgis = mgi.split(",")
			for m in mgis:
				g_anno[m].append(word[3])
	return g_anno


def read_pheno_anno(ifile):
	p_anno = {}
	with open(ifile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) == 3:
				if word[2] != "OBSOLETE.":
					p_anno[word[0]] = word[1]
			elif len(word) == 2:
				#if word[2] != "OBSOLETE.":
				p_anno[word[0]] = word[1]
	return p_anno

def combine_everything(data, gene_dict, g_anno, p_anno, outfile, header):
	output = open(outfile, "w")
	output.write("{}\tPhenotype annotation\n".format(header)),
	for gene in data:
		output.write("{}\t".format(data[gene])),
		if gene in gene_dict:
			if gene_dict[gene] in g_anno: #Now have all pheno_ids
				c = 0
				for pheno_id in g_anno[gene_dict[gene]]:	
					if pheno_id in p_anno:
						if c == 0:
							output.write("{}".format(p_anno[pheno_id])),
						else:
							output.write(",{}".format(p_anno[pheno_id])),
						#m = re.search("lethality", p_anno[pheno_id])
						#m2 = re.search("postnatal", p_anno[pheno_id])
					else:
						if c == 0:
							output.write("No phenotype for gene"),
					c += 1
			else:
				output.write("Not found in phenotype database"),
		else:
			output.write("No MGI ID for gene"),
		output.write("\n"),

def main():
	parser = argparse.ArgumentParser(description='Annotation of gene lists to phenotypic descriptors\n')
	parser.add_argument('-i','--input', help='Input file in tab delimiated format containing ensembl IDs on the first column. Assumes a header is present.', required=True)
	parser.add_argument('-o','--outfile', help='Output file', required=True)
	args = vars(parser.parse_args())

	gene_ids = pkg_resources.resource_filename('pygprofile', 'data/MGI_PhenoGenoMP.rpt')
	pheno_ids = pkg_resources.resource_filename('pygprofile', 'data/VOC_MammalianPhenotype.rpt')

	gene_list, header = read_input(args["input"])
	converted_genes = annotate_ensembl(gene_list)
	gene_anno = read_gene_anno(gene_ids)
	pheno_anno = read_pheno_anno(pheno_ids)
	combine_everything(gene_list, converted_genes, gene_anno, pheno_anno, args["outfile"], header)

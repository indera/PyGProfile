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
import pygenerich
import rpy2.robjects.numpy2ri as rpyn
import numpy as np

def annotate_ensembl(dict_obj):
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	genome="mmusculus_gene_ensembl"
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	values = []
	for key1 in dict_obj:
		values.append(key1)
	C1BM = ro.r.getBM(attributes=StrVector(["ensembl_gene_id", "entrezgene"]), filters="ensembl_gene_id", values=values, mart=ensembl)
	gene = list(C1BM.rx(True,1))
	ent = list(C1BM.rx(True,2))
	data = {}
	for index, g in enumerate(gene):
		data[g] = ent[index]
	return data
	
def spia(gene_vals, pval):
	importr("SPIA")
	result = annotate_ensembl(gene_vals)
	sigs = {}
	univ = []
	for key in result:
		if float(gene_vals[key][1]) < float(pval):
			entrez = str(result[key])
			sigs[entrez] = float(gene_vals[key][0]) #Contains LFC!
		univ.append(result[key])
	sig_genes = _dict_to_namedvector(sigs)
	spia_result = ro.r.spia(de=sig_genes, all=univ, organism="mmu", plots=False)
	output = open("SPIA_result.txt", "w")
	output.write("Name\tID\tPathway Size\tNo. of Diff. Exp. Genes (NDE)\tprobability of NDE\tPertubation of Pathway(tA)\tProbability of tA\tCombination of pNDE and ptA\tFDR\tBonferroni\tStatus\n"),
	for i, go_id in enumerate(spia_result[0]):
		output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(spia_result[0][i], spia_result[1][i],spia_result[2][i], spia_result[3][i],spia_result[4][i],spia_result[5][i],
			spia_result[6][i],spia_result[7][i],spia_result[8][i],spia_result[9][i],spia_result[10][i])),
	output.close()
	ro.r('pdf("SPIA.pdf")')
	ro.r.plotP(spia_result, threshold=0.05)
	ro.r('dev.off()')

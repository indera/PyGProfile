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

def go_stats(gene_vals, pval, outname):
	importr("org.Mm.eg.db")
	importr("GO.db")	
	importr("GOstats")
	result = annotate_ensembl(gene_vals)
	#print result
	sigs = []
	univ = []
	for key in result:
		if gene_vals[key] < float(pval):
			sigs.append(result[key])
		univ.append(result[key])
	bp = ro.r.new('GOHyperGParams', geneIds=sigs, universeGeneIds=univ, ontology='BP', pvalueCutoff=0.001, conditional=False, testDirection='over', annotation="org.Mm.eg.db")
	mf = ro.r.new('GOHyperGParams', geneIds=sigs, universeGeneIds=univ, ontology='MF', pvalueCutoff=0.001, conditional=False, testDirection='over', annotation="org.Mm.eg.db")
	cc = ro.r.new('GOHyperGParams', geneIds=sigs, universeGeneIds=univ, ontology='CC', pvalueCutoff=0.001, conditional=False, testDirection='over', annotation="org.Mm.eg.db")
	hgOver1 = ro.r.hyperGTest(bp)
	hgOver2 = ro.r.hyperGTest(mf)
	hgOver3 = ro.r.hyperGTest(cc)
	result1 = ro.r.summary(hgOver1)
	result2 = ro.r.summary(hgOver2)
	result3 = ro.r.summary(hgOver3)
	output = open("{}_GOstats_BP.txt".format(outname), "w")
	output.write("GOBPID\tPvalue\tOddsRatio\tExpCount\tCount\tSize\tTerm\n"),
	for i, go_id in enumerate(result1[0]):
		output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(result1[0][i], result1[1][i],result1[2][i], result1[3][i],result1[4][i],result1[5][i], result1[6][i])),
	output.close()
	output = open("{}_GOstats_MF.txt".format(outname), "w")
	output.write("GOBPID\tPvalue\tOddsRatio\tExpCount\tCount\tSize\tTerm\n"),
	for i, go_id in enumerate(result2[0]):
		output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(result2[0][i], result2[1][i],result2[2][i], result2[3][i],result2[4][i],result2[5][i], result2[6][i])),
	output.close()
	output = open("{}_GOstats_CC.txt".format(outname), "w")
	output.write("GOBPID\tPvalue\tOddsRatio\tExpCount\tCount\tSize\tTerm\n"),
	for i, go_id in enumerate(result3[0]):
		output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(result3[0][i], result3[1][i],result3[2][i], result3[3][i],result3[4][i],result3[5][i], result3[6][i])),
	output.close()

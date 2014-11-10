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

def custom_pathways(gene_vals, kegg_file, pval):
	importr("KEGGREST")
	importr("org.Mm.eg.db")
	importr("GSEABase")
	result = annotate_ensembl(gene_vals)
	sigs = []
	univ = []
	for key in result:
		if float(gene_vals[key]) < float(pval):
			sigs.append(result[key])
		univ.append(result[key])
	ro.globalenv["sigs"] = sigs
	ro.globalenv["univ"] = univ
	sets = ro.r.getGmt(kegg_file)
	ro.globalenv["sets"] = sets
	ro.r('genes_pathway <- lapply(sets, geneIds)')
	ro.r('names(genes_pathway) <- names(sets)')
	ro.r('hyperg <- Category:::.doHyperGInternal')
	ro.r('''hyperg_test <- function(pathway_genes, significant_genes, all_genes, over=TRUE) {
			white_balls_drawn <- length(intersect(significant_genes, pathway_genes))
			white_balls_in_urn <- length(pathway_genes)
			total_balls_in_urn <- length(all_genes)
			black_balls_in_urn <- total_balls_in_urn - white_balls_in_urn
			balls_pulled <- length(significant_genes)
			hyperg(white_balls_in_urn, black_balls_in_urn, balls_pulled, white_balls_drawn, over) } ''')
	ro.r('pVals_pathway <- t(sapply(genes_pathway, hyperg_test, sigs, univ))')
	ro.r('pVals_pathway <- cbind(rownames(pVals_pathway), pVals_pathway)')
	pvals = ro.r('pVals_pathway')
	vector1=rpyn.ri2numpy(pvals.rx(True,1))
	vector2=rpyn.ri2numpy(pvals.rx(True,2))
	vector3=rpyn.ri2numpy(pvals.rx(True,3))
	vector4=rpyn.ri2numpy(pvals.rx(True,4))
	output = open("Hypergeo_pathways.txt", "w")
	output.write("Pathway\tP-value\tOddsRatio\tExpected\n"),
	for i, j in enumerate(vector1[0]):
		output.write("{}\t{}\t{}\t{}\n".format(j, vector2[0][i],vector3[0][i],vector4[0][i])),
	output.close()

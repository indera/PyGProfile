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

def goprofiles(gene_vals, pval, outname):
	importr("org.Mm.eg.db")
	importr("goProfiles")
	values = []
	for key1 in gene_vals.keys():
		if gene_vals[key1] <= float(pval):
			values.append(key1)
	entrez = annotate_ensembl(values)
	entrezlist = []
	for key1 in entrez.keys():
		entrezlist.append(entrez[key1])
	MF = ro.r.basicProfile(entrezlist, onto="MF", level=2, orgPackage="org.Mm.eg.db")
	BP = ro.r.basicProfile(entrezlist, onto="BP", level=2, orgPackage="org.Mm.eg.db")
	CC = ro.r.basicProfile(entrezlist, onto="CC", level=2, orgPackage="org.Mm.eg.db")
	ro.r('pdf("{}_goProfile.pdf")'.format(outname))
	ro.r('par(mfrow=c(2,2))')
	ro.r.plotProfiles(MF, aTitle="Molecular Function goProfile")
	ro.r.plotProfiles(CC, aTitle="Cellular Component goProfile")
	ro.r.plotProfiles(BP, aTitle="Biological Process goProfile")
	ro.r('dev.off()')
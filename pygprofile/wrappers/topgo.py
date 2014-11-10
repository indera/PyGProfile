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

def run_topgo(gene_vals, pval, outname):
	importr("topGO")
	importr("org.Mm.eg.db")
	importr("biomaRt")
	go_pval = 0.2
	ro.r('''topDiffGenes = function(allScore) { return (allScore < %s)}''' % float(pval)) #Formula for deciding which genes are important
	go_term_type = ["MF", "BP", "CC"]
	topgo_method = "classic" # choice of classic, elim, weight
	for go_type in go_term_type:
		params = {"ontology" : go_type,
			"annot" : ro.r["annFUN.org"],
			"geneSelectionFun" : ro.r["topDiffGenes"],
			"allGenes" : _dict_to_namedvector(gene_vals),
			"mapping" : "org.Mm.eg.db", 
			"ID" : "Ensembl"}
		go_data = ro.r.new("topGOdata", **params)
		results = ro.r.runTest(go_data, algorithm=topgo_method, statistic="fisher")
		scores = ro.r.score(results)
		ro.globalenv["scores"] = scores
		score_names = ro.r('names(scores)')
		num_summarize = min(100, len(score_names))
		results_table = ro.r.GenTable(go_data, elimFisher=results,orderBy="elimFisher", topNodes=num_summarize)
		output = open("{}_topGO_{}.txt".format(outname, go_type), "w")
		output.write("GO.ID\tTerm\tAnnotated\tSignificant\tExpected\telimFisher\n"),
		for index, go_id in enumerate(results_table[0]):
			output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(results_table[0][index], results_table[1][index],results_table[2][index], results_table[3][index],results_table[4][index],results_table[5][index])),
		output.close()
	#Not used but for getting genes to GO
	#goID <- allRes[10, "GO.ID"]
	#gt <- printGenes(GOdata, whichTerms = goID, chip = affyLib, numChar = 40)
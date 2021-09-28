#!/usr/bin/env python

"""
fusion_report_referrals.py

Author: Seemu Ali 
Created: 10.06.2020
Version: 0.0.2

"""

import numpy as np
import pandas as pd 
import sys
import os
import fnmatch

sample_dir=sys.argv[1]
panel_dir=sys.argv[2]
sampleId=sys.argv[3]

tool_output =[ "StarFusion", "Arriba", "ArribaDiscarded", ]
cases = ["IL1", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL8", "IL9", "IL10", "IL11", "IL12", "IL13", "IL14", "IL15", "IL16", "IL17", "IL18", "IL19", "IL20", "IL21", "IL22", "IL23", "IL24", "IL25", "IL26", "IL27", "IL28", "IL29", "IL30"] 


#	*Referrals* 


cases_df=pd.read_csv(panel_dir+'/HaemCases.csv', sep=",", index_col=False)

#def get_genes_list(df, sampleId):
    #case_genes = cases_df[sampleId].to_list()
    #return case_genes 
#get_genes_list(cases_df, sampleId) 



#	*Format Starfusion Report*

fusion_report=pd.read_csv(sample_dir+"/STAR-Fusion/"+sampleId+"_fusionReport.txt", sep="\t", index_col=False)
if (fusion_report.shape[0]!=0):
    genes=fusion_report["Fusion_Name"].str.split("--", expand=True)
    fusion_report["gene1"]=genes[0]
    fusion_report["gene2"]=genes[1]
    fusion_report=pd.DataFrame(fusion_report)
elif (fusion_report.shape[0]==0):
    fusion_report["gene1"]="" 
    fusion_report["gene2"]=""
    fusion_report=pd.DataFrame(fusion_report)


#	*Format Arriba Report*

arriba_report=pd.read_csv(sample_dir+"/Arriba/"+sampleId+"_arriba_fusions.tsv", sep="\t")
arriba_report=arriba_report.rename(columns={'#gene1': 'gene1'})

arriba_discarded=pd.read_csv(sample_dir+"/Arriba/"+sampleId+"_fusions_arriba_discarded.tsv", sep="\t")
arriba_discarded=arriba_discarded.rename(columns={'#gene1': 'gene1'})


#	*Report Generator* 


def report_maker(tool, sampleId, case, genes, results):

	case_dict = { i : case for i in genes}

	results["gene1_mark"] = results["gene1"].map(case_dict) 
	results["gene2_mark"] = results["gene2"].map(case_dict)  
	final_report = results	

	final_report["1-hit"] = np.where((final_report["gene1_mark"] == case) | (final_report["gene2_mark"] == case), 1, 0)  
	final_report["2-hit"] = np.where((final_report["gene1_mark"] == case) & (final_report["gene2_mark"] == case), 1, 0) 
	
	full_report = final_report.loc[final_report["1-hit"] == 1 ]

	if (full_report.shape[0]!=0):
		full_report.to_csv(f"./Results/Full_Reports/{sampleId}_{tool}_Full_Report_unmasked.csv", sep=',')
	if (full_report.shape[0] == 0):
		with open(f"./Results/Full_Reports/Nofusions_unmasked.txt", 'a+') as file:
			file.seek(0)
			data = file.read(100)
			if len(data) > 0 :
        			file.write("\n")
			file.write( f' \n *{tool}: Case-relevant fusions not detected*')
			file.close()
    

for case in cases:
	if case == sampleId:
		genes = cases_df[case].to_list()
		for tool in tool_output:
			if tool == "StarFusion":
				results = fusion_report
				report_maker(tool, sampleId, case, genes, results)
			elif tool == "Arriba":
				results =  arriba_report
				report_maker(tool, sampleId, case, genes, results)
			elif tool == "ArribaDiscarded":
				results = arriba_discarded
				report_maker(tool, sampleId, case, genes, results)

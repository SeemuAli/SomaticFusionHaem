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

sample_dir=sys.argv[1]
panel_dir=sys.argv[2]
sampleId=sys.argv[3]

tool_output =[ "StarFusion", "Arriba", "ArribaDiscarded", ]
referrals = ["Lymphoma", "AML", "CML", "MPN", "ALL", "Myeloma", ] 


#	*Referrals* 


referral_df=pd.read_csv(panel_dir+'/HaemReferrals.csv', sep=",", index_col=False)

lymphoma = referral_df['Lymphoma'].to_list()
AML = referral_df['AML'].to_list()
CML = referral_df['CML'].tolist()
MPN = referral_df['MPN'].tolist()
ALL = referral_df['ALL'].tolist()
myeloma = referral_df['Myeloma'].tolist()


#	*Format Starfusion Report*

fusion_report=pd.read_csv(sample_dir+"/STAR-Fusion/fusionReport/"+sampleId+"_fusionReport.txt", sep="\t", index_col=False)
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


def report_maker(tool, sampleId, referral, genes, results):

	referral_dict = { i : referral for i in genes}

	results["gene1_mark"] = results["gene1"].map(referral_dict) 
	results["gene2_mark"] = results["gene2"].map(referral_dict)  
	final_report = results	

	final_report["1-hit"] = np.where((final_report["gene1_mark"] == referral) | (final_report["gene2_mark"] == referral), 1, 0)  
	final_report["2-hit"] = np.where((final_report["gene1_mark"] == referral) & (final_report["gene2_mark"] == referral), 1, 0) 	
	full_report = final_report.loc[final_report["1-hit"] == 1 ]
	twohit_report = final_report.loc[final_report["2-hit"] == 1 ]
	twohit_report.to_csv(f"./Results/{referral}/{sampleId}_{referral}_{tool}_2Hit_Report.csv", sep=',')
	if (full_report.shape[0]!=0):
		full_report.to_csv(f"./Results/{referral}/Full_Reports/{sampleId}_{referral}_{tool}_Full_Report.csv", sep=',')
	if (twohit_report.shape[0] == 0):
		os.remove(f"./Results/{referral}/{sampleId}_{referral}_{tool}_2Hit_Report.csv")
	if (full_report.shape[0] == 0):
		with open(f"./Results/{referral}/Full_Reports/Nofusions.txt", 'a+') as file:
			file.seek(0)
			data = file.read(100)
			if len(data) > 0 :
        			file.write("\n")
			file.write( f' \n *{tool}: {referral} fusions not detected*')
			file.close()
    

for referral in referrals:
	genes = referral_df[referral].to_list()

	for tool in tool_output:
		if tool == "StarFusion":
			results = fusion_report
			report_maker(tool, sampleId, referral, genes, results)
		elif tool == "Arriba":
			results =  arriba_report
			report_maker(tool, sampleId, referral, genes, results)
		elif tool == "ArribaDiscarded":
			results = arriba_discarded
			report_maker(tool, sampleId, referral, genes, results)





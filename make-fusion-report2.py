#!/usr/bin/env python
# coding: utf-8

"""
make-fusion-report.py

reformats / aggregates STAR-Fusion output \
generates a single text report and and IGV batch file for each sample

Aurthor: Christopher Medway
Created: 12th September 2019
Version: 0.0.1
"""

import os
import argparse


def get_args():
    """
    parses command line arguments
    """
    parser = argparse.ArgumentParser(
        description='generates a single report from STAR-fusion output files'
    )

    parser.add_argument('--sampleId', help='unique identification of sample', required=True)
    parser.add_argument('--seqId', help='run directory ID', required=True)
    parser.add_argument('--panel', help='name of panel', required=True)
    parser.add_argument('--ip', help='required for igv batchfile', required=True)

    args = parser.parse_args()
    return (args)


def make_igv_batchfile(args):
    """
    generates a batch file that can be uploaded into IGV
    :param args:
    :return:
    """
    outputfile = "./fusionReport/" + args.sampleId + "_igv_report.batch"

    path_for_igv = '//' + \
                   args.ip + \
                   '/results/' + \
                   args.seqId + \
                   '/' + \
                   args.panel + \
                   '/' + \
                   args.sampleId + \
                   '/STAR-Fusion/FusionInspector-validate/'

    print(path_for_igv)

    try:
        os.remove(outputfile)
    except OSError:
        pass

    out = open(outputfile, 'a')
    out.write(
        "new" + "\n" +
        "genome " + path_for_igv + "finspector.fa" + "\n" +
        "load " + path_for_igv + "cytoBand.txt" + "\n" +
        "load " + path_for_igv + "finspector.gtf" + "\n" +
        "load " + path_for_igv + "finspector.junction_reads.bam" + "\n" +
        "load " + path_for_igv + "finspector.spanning_reads.bam"
    )

    out.close()


def fusion_report(args):
    """
    compiles a text report from STAR-Fusion input files
    :param args:
    :return:
    """
    outputfile = "./fusionReport/" + args.sampleId + "_fusionReport.txt"
    star_fusion_results_path = "./STAR-Fusion/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv.annotated.coding_effect"

    try:
        os.remove(outputfile)
    except OSError:
        pass

    out = open(outputfile, 'a')
    out.write(
        "Fusion_Name" + "\t" +
        "Split_Read_Count" + "\t" +
        "Spanning_Read_Count" + "\t" +
        "Left_Breakpoint" + "\t" +
        "Right_Breakpoint" + "\t" +
        "CDS_Left_ID" + "\t" +
        "CDS_Left_Range" + "\t" +
        "CDS_Right_ID" + "\t" +
        "CDS_Right_Range" + "\t" +
        "Prot_Fusion_Type" + "\t" +
        "Num_WT_Fragments_Left" + "\t" +
        "Num_WT_Fragments_Right" + "\t" +
        "Fusion_Allelic_Fraction" + "\n")

    with open(star_fusion_results_path) as sf:
        header_line = sf.readline()
        for line in sf:
            sfln = line.split('\t')

            large_anchor_support = sfln[8]
            junc_read_ct = int(sfln[1])
            span_frag_ct = int(sfln[2])

            if large_anchor_support != 'NO_LDAS' and junc_read_ct > 0 and span_frag_ct > 0 and (
                    junc_read_ct + span_frag_ct) > 2:
                fusion = sfln[0]
                left_breakpoint = sfln[5]
                right_breakpoint = sfln[7]
                num_counter_fusion_left = int(sfln[11])
                num_counter_fusion_right = int(sfln[12])

                # caculate the fusion alleleic ratio (FAR)
                fafL = (junc_read_ct + span_frag_ct) / (num_counter_fusion_left + junc_read_ct + span_frag_ct)
                fafR = (junc_read_ct + span_frag_ct) / (num_counter_fusion_right + junc_read_ct + span_frag_ct)
                faf = ((fafL + fafR) / 2)

                out.write(
                    sfln[0] + "\t" +
                    sfln[1] + "\t" +
                    sfln[2] + "\t" +
                    sfln[5] + "\t" +
                    sfln[8] + "\t" +
                    sfln[21] + "\t" +
                    sfln[22] + "\t" +
                    sfln[23] + "\t" +
                    sfln[24] + "\t" +
                    sfln[25] + "\t" +
                    str(num_counter_fusion_left) + "\t" +
                    str(num_counter_fusion_right) + "\t" +
                    str(round(faf, 2)) + "\n"
                )

        out.close()


if __name__ == '__main__':
    os.mkdir("./fusionReport")

    args = get_args()
    fusion_report(args)
    make_igv_batchfile(args)

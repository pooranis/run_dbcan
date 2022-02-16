#!/usr/bin/env python3

import sys
from dbcan.utils import printmsg
import os


_TOOLNAMES = ['diamond', 'hmmer', 'eCAMI']

def make_overview(outDir,prefix,tools=None,use_signalP=False):
    """
    Args:
        outDir: directory containing run_dbcan outputs
        prefix: prefix for run_dbcan files
        tools (dict): dictionary with keys ['diamond', 'hmmer', 'eCAMI'] and values
            bool True or False for which tool output to include in overview. Default
            True for all.
        use_signalP (bool): include signalP outputs (must exist if True or will through error)
    """
    # start Overview
    if not tools:
        tools = dict.fromkeys(_TOOLNAMES, True)
    printmsg(f"Preparing overview table from {', '.join(list(filter(tools.get, tools)))} output...")
    workdir= os.path.join(outDir,prefix)
    # a function to remove duplicates from lists while keeping original order
    def unique(seq):
        exists = set()
        return [x for x in seq if not (x in exists or exists.add(x))]

    # check if files exist. if so, read files and get the gene numbers
    diamond_genes =[]
    hmmer_genes   = []
    eCAMI_genes  =[]
    if tools["diamond"]:
        arr_diamond = open(workdir+"diamond.out").readlines()
        diamond_genes = [arr_diamond[i].split()[0] for i in range(1, len(arr_diamond))] # or diamond_genes = []

    if tools["hmmer"]:
        arr_hmmer = open(workdir+"hmmer.out").readlines()
        hmmer_genes = [arr_hmmer[i].split()[2] for i in range(1, len(arr_hmmer))] # or hmmer_genes = []

    if tools["eCAMI"]:
        arr_eCAMI = open(workdir+"eCAMI.out").readlines()
        eCAMI_genes = [arr_eCAMI[i].split()[0] for i in range(1, len(arr_eCAMI))]# or eCAMI_genes = []

    if use_signalP and (os.path.exists(workdir + "signalp.out")):
        arr_sigp = open(workdir+"signalp.out").readlines()
        sigp_genes = {}
        for i in range (0,len(arr_sigp)):
            row = arr_sigp[i].split()
            sigp_genes[row[0]] = row[4] #previous one is row[2], use Y-score instead from suggestion of Dongyao Li

    ##Catie Ausland edits BEGIN, Le add variable exists or not, remove duplicates from input lists


    if len(eCAMI_genes) > 0:
        if (eCAMI_genes[-1] == None):
            #print('I am in &&&&&&&&&&&&&&&&&&&&&&')
            eCAMI_genes.pop()
            eCAMI_genes = unique(eCAMI_genes)
            if 'hmmer_genes' in locals():
                hmmer_genes.pop()
                hmmer_genes = unique(hmmer_genes)
            if 'diamond_genes' in locals():
                diamond_genes.pop()
                diamond_genes = unique(diamond_genes)
    ## Catie edits END, Le add variable exists or not, remove duplicates from input lists

    # parse input, stroe needed variables
    if tools["diamond"] and (len(arr_diamond) > 1):
        diamond_fams = {}
        for i in range (1,len(arr_diamond)):
            row = arr_diamond[i].split("\t")
            fam = row[1].strip("|").split("|")
            diamond_fams[row[0]] = fam[1:]

    if tools["hmmer"] and (len(arr_hmmer) > 1):
        hmmer_fams = {}
        for i in range (1, len(arr_hmmer)):
            row = arr_hmmer[i].split("\t")
            fam = row[0].split(".")
            fam = fam[0]+"("+row[7]+"-"+row[8]+")"
            if(row[2] not in hmmer_fams):
                hmmer_fams[row[2]] = []
            hmmer_fams[row[2]].append(fam)


    if tools["eCAMI"] and (len(arr_eCAMI) > 1) :
        eCAMI_fams = {}
        for i in range (1,len(arr_eCAMI)):
            row_ori = arr_eCAMI[i].split("\t")
            subfam_names = row_ori[2].split('|')
            fam = row_ori[1].split(':')
            if ' ' in row_ori[0]:
                fams_ID = row_ori[0].split(' ')[0]
            else:
                fams_ID = row_ori[0]

            diam_name = []
            for name in subfam_names:
                if '.' in name:
                    diam_name.append(name.split(":")[0])

            if(fams_ID not in eCAMI_fams):
                eCAMI_fams[fams_ID] = {}
                eCAMI_fams[fams_ID]["fam_name"] = []
                eCAMI_fams[fams_ID]["ec_num"] = []

            eCAMI_fams[fams_ID]["fam_name"].append(fam[0])
            eCAMI_fams[fams_ID]["ec_num"] = diam_name

    #overall table

    all_genes = unique(hmmer_genes+eCAMI_genes+diamond_genes)

    with open(workdir+"overview.txt", 'w+') as fp:
        if use_signalP:
            fp.write("Gene ID\tEC#\tHMMER\teCAMI\tDIAMOND\tSignalp\t#ofTools\n")
        else:
            fp.write("Gene ID\tEC#\tHMMER\teCAMI\tDIAMOND\t#ofTools\n")
        for gene in all_genes:
            csv=[gene]
            num_tools = 0

            if tools["eCAMI"] and arr_eCAMI != None and (gene in eCAMI_genes):
                if eCAMI_fams[gene]["ec_num"] == []:
                    csv.append("-")
                else:
                    csv.append("|".join(eCAMI_fams[gene]["ec_num"]))
            else:
                csv.append("-")

            if tools["hmmer"] and arr_hmmer != None and (gene in hmmer_genes):
                num_tools += 1
                csv.append("+".join(hmmer_fams[gene]))
            else:
                csv.append("-")

            if tools["eCAMI"] and arr_eCAMI != None and (gene in eCAMI_genes):
                num_tools += 1
                csv.append("+".join(eCAMI_fams[gene]["fam_name"]))
            else:
                csv.append("-")

            if tools["diamond"] and arr_diamond != None and (gene in diamond_genes):
                num_tools += 1
                csv.append("+".join(diamond_fams[gene]))
            else:
                csv.append("-")
            if use_signalP:
                if (gene in sigp_genes):
                    csv.append("Y(1-"+sigp_genes[gene]+")")
                else:
                    csv.append("N")
            csv.append(str(num_tools))
            temp = "\t".join(csv) + "\n"
            fp.write(temp)
    printmsg("overview table complete. Saved as "+workdir+"overview.txt")
    # End overview



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='combine tool outputs into overview table')
    req = parser.add_argument_group('required arguments')
    req.add_argument('-d', '--outDir', help='Directory containing result of run_dbcan', required=True)
    parser.add_argument('-p', '--prefix', help='Input files prefix from result of run_dbcan', default='')
    parser.add_argument('--tools', '-t', action='append', choices=['hmmer', 'diamond', 'eCAMI', 'all'], help='Choose a combination of tools to run. Can specify 2 by repeating option, e.g. `-t hmmer -t eCAMI`. (default: all)')
    parser.add_argument('--use_signalP', action='store_true', help='Include signalP output.')
    parser._action_groups.reverse()
    args = parser.parse_args()

    ## make tools dictionary
    tools = dict.fromkeys(_TOOLNAMES, False)
    if args.tools is None or 'all' in args.tools:
        args.tools = _TOOLNAMES

    tools.update(dict.fromkeys(args.tools, True))

    make_overview(outDir = args.outDir ,prefix = args.prefix,tools=tools,use_signalP=False)

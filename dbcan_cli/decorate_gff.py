#!/usr/bin/env python3
"""
Decorates input gff file with annotations from run_dbcan.
For usage::

    make_gff.py --help
    
"""


import os
from importlib.metadata import version
import argparse
import errno
from dbcan.utils import printmsg

CONST_SOURCE = "run_dbcan_"+version('dbcan')

def write_gff(auxFile, outputgff, source = CONST_SOURCE,
             cazyme_genes = None, tf_genes = None, tp_genes = None, stp_genes = None):
    """write gff file

    - auxFile: input GFF file
    - outputgff: output GFF file
    - source: what to put in source field of gff
    - cazyme_genes, tf_genes, tp_genes, stp_genes: dictionaries mapping from ID in gff to whatever annotation

    """
    (cazyme_genes, tf_genes, tp_genes, stp_genes) = [x or {} for x in [cazyme_genes, tf_genes, tp_genes, stp_genes]]
    printmsg(f"Writing output gff file {outputgff}")
    with open(auxFile) as f:
        with open(outputgff, 'w') as out:
            for line in f:
                if not line.startswith("#"):
                    row = line.rstrip().split('\t')
                    if row[2] == "CDS":
                        note = row[8].strip().rstrip(";").split(";")
                        ID = ""
                        notes = {}
                        for x in note:
                            temp = x.split('=')
                            notes[temp[0]] = temp[1]

                        if "ID" in notes:
                            ID = notes["ID"]
                        elif "Name" in notes:
                            ID = notes["Name"]
                        else:
                            continue

                        parent = ID
                        if ID in cazyme_genes:
                            row[2] = "CAZyme"
                            gffdb = ";DB="+'|'.join(cazyme_genes[ID])
                            ID += "_"+"cazyme"
                        elif ID in tf_genes:
                            row[2] = "TF"
                            gffdb = ";DB="+tf_genes[ID]
                            ID += "_"+"tf"
                        elif ID in tp_genes:
                            row[2] = "TC"
                            gffdb = ";DB="+tp_genes[ID]
                            ID += "_"+"tc"
                        elif ID in stp_genes:
                            row[2] = "STP"
                            gffdb = ";DB=" + stp_genes[ID]
                            ID += "_"+"stp"
                        else:
                            continue
                        row[8] = "ID=" + ID + gffdb
                        row[1] = source
                        if "Name" in notes:
                            row[8] += ";Name="+notes["Name"]
                        row[8] += ';Parent=' + parent + ";locus_tag=" + parent
                        out.write(line)
                        out.write('\t'.join(row)+'\n')


def get_cgc_genes(out_dir, out_pre):
    """get genes from CGC finder step output

    Args:
        out_dir (str): directory containing output of CGC finder
        out_pre (str): prefix of filenames (can be '')
    
    Returns:
        tuple of 3 dictionaries mapping from sequence IDs to output of search with
        particular DB annotation info.
        In order:
        tf_genes (dict): tf-1 and tf-2 hmm dbs
        tp_genes (dict): (sic) tcdb diamond db
        stp_genes (dict): stp hmm db
    """
    outPath = os.path.join(out_dir, out_pre)
    tf_genes = {}
    tp_genes = {}
    stp_genes = {}

    try:
        with open(f"{outPath}tf-1.out") as f:
            printmsg(f"getting TF genes from {outPath}tf-1.out")
            for line in f:
                row = line.rstrip().split('\t')
                row[0] = "DBD-Pfam|" + row[0]
                if not row[2] in tf_genes:
                    tf_genes[row[2]] = row[0]
                else:
                    tf_genes[row[2]] += ',' + row[0]

        with open(f"{outPath}tf-2.out") as f:
            printmsg(f"getting TF genes from {outPath}tf-2.out")
            for line in f:
                row = line.rstrip().split('\t')
                row[0] = "DBD-SUPERFAMILY|" + row[0]
                if not row[2] in tf_genes:
                    tf_genes[row[2]] = row[0]
                else:
                    tf_genes[row[2]] += ',' + row[0]

        with open(f"{outPath}tp.out") as f:
            printmsg(f"getting TP/TC genes from {outPath}tp.out")
            for line in f:
                row = line.rstrip().split('\t')
                if not row[0] in tp_genes:
                    tp_genes[row[0]] = row[1]
                else:
                    tp_genes[row[0]] += ','+row[1]

        with open(f"{outPath}stp.out") as f:
            printmsg(f"getting STP genes from {outPath}stp.out")
            for line in f:
                row = line.rstrip().split('\t')
                row[0] = "STP|" + row[0]
                if not row[2] in stp_genes:
                    stp_genes[row[2]] = row[0]
                else:
                    stp_genes[row[2]] += ',' + row[0]
    except FileNotFoundError as e:
        e.strerror = "Cannot find output of run_dbcan cgc finder step. " + e.strerror
        raise e

    return(tf_genes, tp_genes, stp_genes)

def get_cazyme_genes(out_dir, out_prefix):
    """get cazyme genes from run_dbcan output"""
    outPath = os.path.join(out_dir, out_prefix)
    tools = {}
    cazyme_genes = {}
    dia = set()
    hmm = set()
    eca = set()
    diamondpath = outPath+'diamond.out'
    if os.path.isfile(diamondpath):
        tools["diamond"] = True
        printmsg("getting cazyme genes from", diamondpath)
        with open(diamondpath) as f:
            next(f)
            for line in f:
                row = line.rstrip().split('\t')
                dia.add(row[0])
                if row[0] not in cazyme_genes:
                    cazyme_genes[row[0]] = set()
                cazyme_genes[row[0]].update(set(row[1].strip("|").split('|')[1:]))
    hmmerpath = outPath+'hmmer.out'
    if os.path.isfile(hmmerpath):
        tools["hmmer"] = True
        printmsg("getting cazyme genes from", hmmerpath)
        with open(hmmerpath) as f:
            next(f)
            for line in f:
                row = line.rstrip().split('\t')
                hmm.add(row[2])
                if row[2] not in cazyme_genes:
                    cazyme_genes[row[2]] = set()
                cazyme_genes[row[2]].add(row[0].split('.hmm')[0])
    ecamipath = outPath+'eCAMI.out'
    if os.path.isfile(ecamipath):
        tools["eCAMI"] = True
        printmsg("getting cazyme genes from", ecamipath)
        with open(ecamipath) as f:
            next(f)
            for line in f:
                row_ori = line.rstrip().split('\t')
                if ' ' in row_ori[0]:
                    fams_ID = row_ori[0].split(' ')[0]
                else:
                    fams_ID = row_ori[0]
                eca.add(fams_ID)
                if fams_ID not in cazyme_genes:
                    cazyme_genes[fams_ID] = set()
                cazyme_genes[fams_ID].add(row_ori[1].split(':')[0])

    if list(tools.values()).count(True) == 0:
        raise FileNotFoundError(errno.ENOENT, f"No cazyme outputs ({diamondpath}, {hmmerpath}, or {ecamipath}) of dbcan can be found.")

    if list(tools.values()).count(True) > 1:
        temp1 = hmm.intersection(eca)
        # print(hmm, 'This intersection  hmm')
        temp2 = hmm.intersection(dia)
        # print(dia, 'This intersection  dia')
        temp3 = dia.intersection(eca)
        # print(eca, 'This intersection  eca')
        cazyme = temp1.union(temp2, temp3)
    else:
        cazyme = hmm.union(dia, eca)

    cazyme_genes = {key: cazyme_genes[key] for key in cazyme}
    return(cazyme_genes)

def make_gff(input_gff, output_gff, in_dir, in_prefix, cgc_genes=False, source=CONST_SOURCE):
    """make gff file

    Args:
        in_dir: Input directory containing result of run_dbcan
        in_prefix: Input files prefix from result of run_dbcan
        cgc_genes:bool: Use result of CGC Finder step.
        other fields see :func:`write_gff`
    """
    cazyme_genes = get_cazyme_genes(in_dir, in_prefix)
    (tf_genes, tp_genes, stp_genes) = [None, None, None]
    if cgc_genes:
        (tf_genes, tp_genes, stp_genes) = get_cgc_genes(in_dir, in_prefix)
    write_gff(auxFile = input_gff, source = source, outputgff = output_gff,
             cazyme_genes = cazyme_genes, tf_genes = tf_genes, tp_genes = tp_genes, stp_genes = stp_genes)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='decorate gff file with run_dbcan outputs')
    req = parser.add_argument_group('required arguments')
    req.add_argument('-i', '--input_gff', help='User input GFF file to be decorated/annotated.', required=True)
    req.add_argument('-o', '--output_gff', help='Output GFF file.', required=True)
    req.add_argument('-d', '--in_dir', help='Input directory containing result of run_dbcan', required=True)
    parser.add_argument('-p', '--in_pre', help='Input files prefix from result of run_dbcan', default='')
    parser.add_argument('-c', '--cgc_genes', help='Use result of CGC Finder step.', action='store_true')
    parser.add_argument('-s', '--source', help='Source field for GFF (default: %(default)s).', default=CONST_SOURCE)
    parser._action_groups.reverse()
    args = parser.parse_args()

    make_gff(args.input_gff, args.output_gff, args.in_dir, args.in_pre, cgc_genes = args.cgc_genes, source = args.source)

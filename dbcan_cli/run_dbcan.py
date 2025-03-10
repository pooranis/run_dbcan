#!/usr/bin/env python3
#########################################################
# dbCAN2 Driver Script (Stand Alone Version)
#
# Written by Tanner Yohe in the Yin Lab at NIU
# Revised by Le Huang in the Zhang Lab at NKU
# Updated by Mohamad Majd Raslan in the Yin Lab at NIU
# Updated by Wei Li created table
# Updated by Le Huang at NKU
# Updated by Qiwei Ge in Dr.Yin's Lab at UNL
# Updated by Poorani Subramanian in BCBB at NIAID
# Last updated 12 Feb 2022
# updated information[Qiwei Ge]: 1. Hotpep has been removed, added eCAMI tool. 2. cgc out reformatting. 3. Fixed issues multiple GT2s.
# Accepts user input
# Predicts genes if needed
# Runs input against HMMER, DIAMOND, and Hotpep
# Optionally predicts CGCs with CGCFinder
# Creats an overview table using output files from core
# tools from Hotpep.out,hmmer.out and diamond.out
##########################################################
import subprocess
import glob
import os
import argparse
import sys
import shutil
from dbcan.utils import simplify_output, cgc_finder, printmsg
from dbcan_cli.hmmscan_parser import hmmscan_parse
import dbcan_cli.decorate_gff as decorate_gff
from dbcan_cli.overview import make_overview
from dbcan.eCAMI import eCAMI_config, eCAMI_main

'''
def some functions
'''

def runHmmScan(outPath, hmm_cpu, dbDir, hmm_eval, hmm_cov, db_name, outputFile=None):
    """run hmmscan with {outPath}uniInput query file against dbDir/db_name
    if hmm_cpu > 1, uniInput files should already be split into hmm_cpu parts!
    if outputFile=None, outputFile=f"{outPath}{db_name.replace('.hmm', '')}.out"
    """
    if (hmm_cpu > 1):
        printmsg("Using {hmm_cpu} processes for hmmscan.")
        faa = glob.glob(f"{outPath}uniInput.split/*.faa")
        ps = []
        for i in range(len(faa)):
            cmd = f"hmmscan --domtblout {outPath}uniInput.split/{i}.h.{db_name}.out --noali --cpu 1 -o /dev/null {os.path.join(dbDir, db_name)}.hmm {faa[i]}"
            printmsg(cmd)
            cmdp = subprocess.Popen(cmd, shell=True)
            ps.append(cmdp)
        cmdoutput = [p.wait() for p in ps]
        printmsg(cmdoutput)
        printmsg(f"doing cat {outPath}h{db_name}.out")
        subprocess.run(f"cat {outPath}uniInput.split/*.h.{db_name}.out >{outPath}h{db_name}.out", shell=True)
        [os.remove(splitfile) for splitfile in glob.glob(f"{outPath}uniInput.split/*.h.{db_name}.out")]
    else:
        subprocess.run(['hmmscan', '--domtblout', '%sh%s.out' % (outPath, db_name), '--cpu', str(hmm_cpu), '-o', '/dev/null', '%s%s' % (dbDir,db_name), '%suniInput' % outPath])
    outputFile = outputFile or f"{outPath}{db_name.replace('.hmm', '')}.out"
    hmmscan_parse(f"{outPath}h{db_name}.out", float(hmm_eval), float(hmm_cov), outputFile=outputFile)
    if os.path.exists(f'{outPath}h{db_name}.out'):
        os.remove(f'{outPath}h{db_name}.out')


def cli_main():

    parser = argparse.ArgumentParser(description='dbCAN2 Driver Script')
    parser.add_argument('inputFile', help='User input file. Must be in FASTA format.')
    parser.add_argument('inputType', choices=['protein', 'prok', 'meta'], #protein=proteome, prok=prokaryote nucleotide, meta=metagenome nucleotide
                        help='Type of sequence input. protein=proteome; prok=prokaryote; meta=metagenome')
    parser.add_argument('--tools', '-t', action='append', choices=['hmmer', 'diamond', 'eCAMI', 'all'], help='Choose a combination of tools to run. Can specify more than one by repeating option. (default: all)')
    ## db
    parser.add_argument('--db_dir', default="db", help='Database directory (default: %(default)s)')
    parser.add_argument('--dbCANFile',default="dbCAN.txt", help='Indicate the file name of HMM database such as dbCAN.txt, please use the newest one from dbCAN2 website. (default: %(default)s)')
    ## diamond
    parser.add_argument('--dia_eval', default=1e-102,type=float, help='DIAMOND E Value (default: %(default)s)')
    parser.add_argument('--dia_cpu', default=1, type=int, help='Number of CPU cores that DIAMOND is allowed to use (default: %(default)s)')
    ## hmmer
    parser.add_argument('--hmm_eval', default=1e-15, type=float, help='HMMER E Value (default: %(default)s)')
    parser.add_argument('--hmm_cov', default=0.35, type=float, help='HMMER Coverage val (default: %(default)s)')
    parser.add_argument('--hmm_cpu', default=1, type=int, help='Number of CPU cores that HMMER is allowed to use for hmmer step or CGC tf.hmm and stp.hmm steps. (default: %(default)s)')
    # eCAMI
    parser.add_argument('--eCAMI_kmer_db', default="CAZyme",type=str, help="Change n_mer directories path for prediction (default: %(default)s)")
    parser.add_argument('--eCAMI_k_mer', default=8, type=int, help="Peptide length for prediction (default: %(default)s)")
    parser.add_argument('--eCAMI_jobs', default=1, type=int, help='Number of processor for use for prediction (default: %(default)s)')
    parser.add_argument('--eCAMI_important_k_mer_number', default=5, type=int, help="Minimum number of n_mer for prediction (default: %(default)s)")
    parser.add_argument('--eCAMI_beta', default=2, type=float, help="Minimum sum of percentage of frequency of n_mer for prediction (default: %(default)s)")
    # cluster hmm
    parser.add_argument('--tf_eval', default=1e-4, type=float, help='tf.hmm HMMER E Value (default: %(default)s)')
    parser.add_argument('--tf_cov', default=0.35, type=float, help='tf.hmm HMMER Coverage val (default: %(default)s)')
##    parser.add_argument('--tf_cpu', default=1, type=int, help='tf.hmm Number of CPU cores that HMMER is allowed to use')
    parser.add_argument('--stp_eval', default=1e-4, type=float, help='stp.hmm HMMER E Value (default: %(default)s)')
    parser.add_argument('--stp_cov', default=0.3, type=float, help='stp.hmm HMMER Coverage val (default: %(default)s)')
##    parser.add_argument('--stp_cpu', default=1, type=int, help='stp.hmm Number of CPU cores that HMMER is allowed to use')
    parser.add_argument('--cluster', '-c', help='Predict CGCs via CGCFinder. This argument requires an auxillary locations file if a protein input is being used')
    parser.add_argument('--cgc_dis', default=2, help='CGCFinder Distance value (default: %(default)s)')
    parser.add_argument('--cgc_sig_genes', default='tp', choices=['tp', 'tf','all'], help='CGCFinder Signature Genes value (default: %(default)s)')
    parser.add_argument('--use_signalP', action='store_true', help='Use signalP or not, remember, you need to setup signalP tool first. Because of signalP license, Docker version does not have signalP.')
    parser.add_argument('--gram', '-g', choices=["p","n","all"], default="all", help="Choose gram+(p) or gram-(n) for proteome/prokaryote nucleotide, which are params of SignalP, only if user use signalP")
    # output
    parser.add_argument('--out_pre', default="", help='Output files prefix')
    parser.add_argument('--out_dir', default="output", help='Output directory - full path works best (default: %(default)s)')
    args = parser.parse_args()

    ####
    #run_dbcan.py [inputFile] [inputType]
    ####

    ##########################
    # Begin Setup and Input Checks

    dbDir = args.db_dir
    prefix = args.out_pre
    outDir = args.out_dir

    if not dbDir.endswith("/") and len(dbDir) > 0:
        dbDir += "/"

    if not outDir.endswith("/") and len(outDir) > 0:
        outDir += "/"

    outPath = outDir + prefix
    auxFile = ""
    inputFile = args.inputFile
    inputType = args.inputType
    find_clusters = False
    if args.cluster != None:
        find_clusters = True
        if inputType == "protein":
            auxFile = args.cluster
        else:
            auxFile = '%sprodigal.gff'%outPath

    if not os.path.isdir(dbDir):
        sys.exit(f"ERROR: The database directory {dbDir} does not exist")

    if not os.path.isfile(os.path.join(dbDir,'CAZy.dmnd')):
        sys.exit("ERROR: No CAZy DIAMOND database found. \
        Please make sure that your CAZy DIAMOND databased is named 'CAZy.dmnd' and is located in your database directory")

    if not os.path.isfile(os.path.join(dbDir, args.dbCANFile)):
        sys.exit(f"ERROR: No dbCAN HMM database found. \
        Please make sure that your dbCAN HMM database is named {args.dbCANFile} or the newest one, has been through hmmpress, and is located in your database directory {dbDir}")

    if not os.path.isdir(outDir):
        subprocess.call(['mkdir', outDir])

    if find_clusters and inputType == "protein":
        if len(auxFile) > 0:
            printmsg(auxFile)
            if not os.path.isfile(auxFile):
                    sys.exit(f"ERROR: It seems that the auxillary filename {auxFile} that you provided does not exist, or is not a file")
        else:
            sys.exit("ERROR: Please provide an auxillary input file with the position of each gene. This file can either be in BED or GFF format")
    tools = [True, True, True] #DIAMOND, HMMER, eCAMI
    args.tools = args.tools or ['all'] #default
    if 'all' not in args.tools:
        if 'diamond' not in args.tools:
            tools[0] = False
        if 'hmmer' not in args.tools:
            tools[1] = False
        if 'eCAMI' not in args.tools:
            tools[2] = False


    # End Setup and Input Checks
    #########################
    #########################
    # Begin Gene Prediction Tools
    if inputType == 'prok':
        subprocess.run(['prodigal', '-i', inputFile, '-a', '%suniInput'%outPath, '-o', '%sprodigal.gff'%outPath, '-f', 'gff', '-q'], check=True)
    if inputType == 'meta':
        subprocess.run(['prodigal', '-i', inputFile, '-a', '%suniInput'%outPath, '-o', '%sprodigal.gff'%outPath, '-f', 'gff', '-p', 'meta','-q'], check=True)
    #Proteome
    if inputType == 'protein':
        subprocess.run(['cp', inputFile, '%suniInput'%outPath], check=True)

    ## check for hmm_cpu and split files if > 1
    if (tools[1] or find_clusters) and (args.hmm_cpu > 1):
        printmsg("splitting input for hmmscan parallel")
        if os.path.isdir(f"{outPath}uniInput.split"):
            printmsg(f"WARNING: overwriting {outPath}uniInput.split")
            shutil.rmtree(f"{outPath}uniInput.split")
        subprocess.run(f"seqkit split2 --force -p {args.hmm_cpu} {outPath}uniInput ", shell=True, check=True)
        fs=glob.glob(f"{outPath}uniInput.split/*uniInput*")
        [os.rename(splitfile, splitfile + ".faa") for splitfile in fs]
        ## check number of split files
        if len(fs) < hmm_cpu:
            hmm_cpu=len(fs)
            printmsg(f"WARNING: hmmscan will only use {hmm_cpu} processes as there are only {hmm_cpu} split files, {fs}. (maybe because the number of seqs is < {hmm_cpu})")
        del fs


    # End Gene Prediction Tools
    #######################
    # Begin SignalP
    if args.use_signalP:
        printmsg("***************************0. SIGNALP start*************************************************\n\n", begin="\n\n")
        if args.gram == "p" or args.gram=="all":
            signalpos = subprocess.Popen('signalp -t gram+ %suniInput > %ssignalp.pos' % (outPath, outPath), shell=True)
        if args.gram == "n" or args.gram == "all":
            signalpneg = subprocess.Popen('signalp -t gram- %suniInput > %ssignalp.neg' % (outPath, outPath), shell=True)

    # End SignalP
    #######################
    # Begin Core Tools

    if tools[0]:
        # diamond blastp -d db/CAZy -e 1e-102 -q output_EscheriaColiK12MG1655/uniInput -k 1 -p 2 -o output_EscheriaColiK12MG1655/diamond1.out -f 6
        printmsg("***************************1. DIAMOND start*************************************************\n\n", begin="\n\n")
        os.system('diamond blastp -d %s -e %s -q %suniInput -k 1 -p %d -o %sdiamond.out -f 6'%(os.path.join(dbDir, "CAZy"), str(args.dia_eval), outPath, args.dia_cpu, outPath))
        # diamond = Popen(['diamond', 'blastp', '-d', '%sCAZy.dmnd' % dbDir, '-e', str(args.dia_eval), '-q', '%suniInput' % outPath, '-k', '1', '-p', str(args.dia_cpu), '-o', '%sdiamond.out'%outPath, '-f', '6'])
        printmsg("***************************1. DIAMOND end***************************************************\n\n", begin="\n\n")


    if tools[1]:
        printmsg("***************************2. HMMER start*************************************************\n\n", begin="\n\n")

        runHmmScan(outPath, args.hmm_cpu, dbDir, args.hmm_eval, args.hmm_cov, args.dbCANFile, outputFile=f"{outPath}hmmer.out")

        printmsg("***************************2. HMMER end***************************************************\n\n", begin="\n\n")

        with open(f"{outPath}hmmer.out", "r+") as f:
            text = f.read()
        os.remove(f"{outPath}hmmer.out")
        text = text.split('\n')
        if '' in text:
            text.remove('')
        with open(f"{outPath}hmmer.out", 'a') as f:
            for i in range(len(text)):
                if 'GT2_' in text[i]:
                    profile = text[i].split('\t')[0].split('.')[0]
                    text[i] = text[i].replace(profile,'GT2')
                    f.write(text[i]+'\n')
        if os.path.exists(f"{outPath}h.out"):
            os.remove(f"{outPath}h.out")

    if tools[2]:
        printmsg("***************************3. eCAMI start***************************************************\n\n", begin="\n\n")
        ecami_config = eCAMI_config(
            db_type=args.eCAMI_kmer_db,
            input = f"{outPath}uniInput",
            output= f"{outPath}eCAMI.out",
            k_mer = args.eCAMI_k_mer,
            jobs = args.eCAMI_jobs,
            important_k_mer_number = args.eCAMI_important_k_mer_number,
            beta = args.eCAMI_beta
        )
        eCAMI_main(ecami_config)
        # os.system('python eCAMI/prediction.py -input %suniInput -kmer_db eCAMI/%s -output %seCAMI.out -k_mer %s -jobs %s -important_k_mer_number %s -beta %s' % (outPath, str(args.eCAMI_kmer_db), outPath, str(args.eCAMI_k_mer), str(args.eCAMI_jobs), str(args.eCAMI_important_k_mer_number),str(args.eCAMI_beta)))
        printmsg("***************************3. eCAMI end***************************************************\n\n", begin="\n\n")
    # End Core Tools
    ########################
    # Begin Adding Column Headers

    if tools[2]:
        with open(outPath+'eCAMI.out') as f:
            with open(outPath+'temp', 'w') as out:
                out.write('protein_name\tfam_name:group_number\tsubfam_name_of_the_group:subfam_name_count\n')
                # for line in f:
                for count, line in enumerate(f):
                    if count % 2 == 0:
                        more_information = line.split(">")
                        out.write(more_information[1])
        subprocess.call(['mv', outPath+'temp', outPath+'eCAMI.out'])

    if tools[1]:
        with open(outDir+prefix+'hmmer.out') as f:
            with open(outDir+prefix+'temp', 'w') as out:
                out.write('HMM Profile\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n')
                for line in f:
                    out.write(line)
        subprocess.call(['mv', outDir+prefix+'temp', outDir+prefix+'hmmer.out'])

    if tools[0]:
        with open(outDir+prefix+'diamond.out') as f:
            with open(outDir+prefix+'temp', 'w') as out:
                out.write('Gene ID\tCAZy ID\t% Identical\tLength\tMismatches\tGap Open\tGene Start\tGene End\tCAZy Start\tCAZy End\tE Value\tBit Score\n')
                for line in f:
                    out.write(line)
        subprocess.call(['mv', outDir+prefix+'temp', outDir+prefix+'diamond.out'])

    # End Adding Column Headers
    ########################
    # Begin CGCFinder

    if find_clusters:
        printmsg("*****************************CGC-Finder start************************************")

    ########################
    # Begin TF,TP, STP prediction
        '''
        previous tf uses diamond, Now tf-1 and tf-2 uses hmmer
        tf hmmer
        '''
        #call(['diamond', 'blastp', '-d', dbDir+'tf_v1/tf.dmnd', '-e', '1e-10', '-q', '%suniInput' % outPath, '-k', '1', '-p', '1', '-o', outDir+prefix+'tf.out', '-f', '6'])
        runHmmScan(outPath, args.hmm_cpu, dbDir, str(args.tf_eval), str(args.tf_cov), "tf-1.hmm")
        runHmmScan(outPath, args.hmm_cpu, dbDir, str(args.tf_eval), str(args.tf_cov), "tf-2.hmm")
        '''
        stp hmmer
        '''
        runHmmScan(outPath, args.hmm_cpu, dbDir, str(args.stp_eval), str(args.stp_cov), "stp.hmm")

        '''
        tp diamond
        '''
        subprocess.call(['diamond', 'blastp', '-d', dbDir+'tcdb.dmnd', '-e', '1e-10', '-q', '%suniInput' % outPath, '-k', '1', '-p', '1', '-o', outPath+'tp.out', '-f', '6'])

        try:
            (tf_genes, tp_genes, stp_genes) = decorate_gff.get_cgc_genes(outDir, prefix)
        except FileNotFoundError as e:
            printmsg(e)
            (tf_genes, tp_genes, stp_genes) = [{}, {}, {}]
        # # End TF and TP prediction
        ##########################
        # Begine CAZyme Extraction

        try:
            cazyme_genes = decorate_gff.get_cazyme_genes(outDir, prefix)
        except FileNotFoundError as e:
            printmsg(e)
            cazyme_genes = {}
        cazyme = set(list(cazyme_genes.keys()))
        # End CAZyme Extraction
        ######################
        # Begin GFF preperation

        if inputType == "prok" or inputType == "meta":   #use Prodigal GFF output
            with open(outDir+prefix+'prodigal.gff') as f:
                with open(outDir+prefix+'cgc.gff', 'w') as out:
                    for line in f:
                        if not line.startswith("#"):
                            row = line.rstrip().rstrip(";").split('\t')
                            num = row[-1].split(";")[0].split('_')[-1]
                            gene = row[0] + '_' + num
                            row[8] = ""
                            if gene in cazyme:
                                row[2] = "CAZyme"
                                row[8] = "DB="+'|'.join(cazyme_genes[gene])
                            elif gene in tf_genes:
                                row[2] = "TF"
                                row[8] = "DB="+tf_genes[gene]
                            elif gene in tp_genes:
                                row[2] = "TC"
                                row[8] = "DB="+tp_genes[gene]
                            elif gene in stp_genes:
                                row[2] = "STP"
                                row[8] = "DB="+stp_genes[gene]
                            row[8] += ";ID="+gene
                            out.write('\t'.join(row)+'\n')
        else:  #user provided GFF/BED file
            gff = False
            with open(auxFile) as f:
                for line in f:
                    if not line.startswith('#'):
                        if len(line.split('\t')) == 9:
                            gff = True
                            break
            if gff:  #user file was in GFF format
                decorate_gff.write_gff(auxFile = auxFile, outputgff = outDir+prefix+'cgc.gff',
                                   cazyme_genes = cazyme_genes, tf_genes = tf_genes,
                                   tp_genes = tp_genes, stp_genes = stp_genes)
            else:  #user file was in BED format
                with open(auxFile) as f:
                    with open(outDir+prefix+'cgc.gff', 'w') as out:
                        for line in f:
                            if line.startswith("track"):
                                continue
                            row = line.rstrip().rstrip(";").split('\t')
                            outrow = ['.','.','.','.','.','.','.','.','']
                            gene = row[1]
                            if gene in cazyme:
                                outrow[2] = 'CAZyme'
                                outrow[8] = "DB="+'|'.join(cazyme_genes[gene])
                            elif gene in tf_genes:
                                outrow[2] = 'TF'
                                outrow[8] =  "DB="+tf_genes[gene]
                            elif gene in tp_genes:
                                outrow[2] = 'TC'
                                outrow[8] = "DB="+tp_genes[gene]
                            elif gene in stp_genes:
                                outrow[2] = 'STP'
                                outrow[8] = "DB=" + stp_genes[gene]
                            else:
                                outrow[2] = 'CDS'
                            outrow[0] = row[0]
                            outrow[3] = row[2]
                            outrow[4] = row[3]
                            outrow[6] = row[4]
                            outrow[8] += ";ID="+gene
                            out.write('\t'.join(outrow)+'\n')
        # End GFF
        ####################
        # Begin CGCFinder call

        # call(['CGCFinder.py', outDir+prefix+'cgc.gff', '-o', outDir+prefix+'cgc.out', '-s', args.cgc_sig_genes, '-d', str(args.cgc_dis)])
        cgc_finder(outDir+prefix+'cgc.gff', args.cgc_dis, args.cgc_sig_genes, outDir+prefix+'cgc.out')
        simplify_output(outDir+prefix+'cgc.out')
        printmsg("**************************************CGC-Finder end***********************************************")
        # End CGCFinder call
        # End CGCFinder
        ####################

    ## cleanup split input files if they exist as hmmer will not be run after this
    if os.path.exists(f"{outPath}uniInput.split"):
        shutil.rmtree(f"{outPath}uniInput.split", ignore_errors=True)

    # Begin SignalP combination
    if args.use_signalP:
        printmsg("Waiting on signalP")
        with open(outDir+prefix+'temp', 'w') as out:
            if args.gram == "all" or args.gram =="p":
                signalpos.wait()
                printmsg("SignalP pos complete")

                with open(outDir+prefix+'signalp.pos') as f:
                    for line in f:
                        if not line.startswith('#'):
                            row = line.split(' ')
                            row = [x for x in row if x != '']
                            if row[9] == 'Y':
                                out.write(line)
                subprocess.call(['rm', outDir+prefix+'signalp.pos'])
            if args.gram == "all" or args.gram == "n":
                signalpneg.wait()
                printmsg("SignalP neg complete")
                with open(outDir+prefix+'signalp.neg') as f:
                    for line in f:
                        if not line.startswith('#'):
                            row = line.split(' ')
                            row = [x for x in row if x != '']
                            if row[9] == 'Y':
                                out.write(line)
                subprocess.call(['rm', outDir+prefix+'signalp.neg'])
        subprocess.call('sort -u '+outDir+prefix+'temp > '+outDir+prefix+'signalp.out', shell=True)
        subprocess.call(['rm', outDir+prefix+'temp'])

    # End SignalP combination
    #######################
    # #######################
    # # start Overview
    make_overview(outDir, prefix, tools = dict(zip(['diamond', 'hmmer', 'eCAMI'], tools)), use_signalP = args.use_signalP)
    # # End overview


if __name__ == '__main__':
    cli_main()

#!/usr/bin/env python3
##########################################################
# hmmscan parser for dbCAN meta server
#
# Based off the hmmscan parser used in the dbCAN server,
# written by Dr. Yin
#
# Written by Tanner Yohe under the supervision
# of Dr. Yin in the YinLab at NIU.
#
# Updated by Le Huang from tips the contributor WATSON Mick <mick.watson@roslin.ed.ac.uk>,
# Thank you!
#
# INPUT
# python hmmscan-parser-dbCANmeta.py [inputFile] [eval] [coverage]
# eval and coverage are optional, inputFile is required
# -updating info:
# -adds pid for every subprocess to make codes robust.
# Last updated: 1/10/19
###########################################################


from subprocess import run
import sys
import os

def hmmscan_parse(inputFile, eval, coverage, outputFile=None):
    """
    Args:
        inputFile (filename): output of ``hmmscan --domtblout``
        eval (float): e-value threshold
        coverage (float): coverage threshold in (0,1]
        outputFile (filename): output filename, if None, will use STDOUT
    """
    tmpfile = "temp." + str(os.getpid())

    run("cat "+inputFile+"  | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_.\"\n\";}}' > " + tmpfile, shell=True)

    if outputFile:
        outf = open(outputFile, 'w')
    else:
        outf = sys.stdout


    with open(tmpfile) as f:
        for line in f:
            row = line.rstrip().split('\t')
            row.append(float(int(row[6])-int(row[5]))/int(row[1]))
            if float(row[4]) <= eval and float(row[-1]) >= coverage:
                print('\t'.join([str(x) for x in row]), file=outf)

    os.remove(tmpfile)
    if outputFile:
        outf.close()


if __name__ == "__main__":
    eval = 1e-15
    coverage = 0.35

    if len(sys.argv) > 1:
        inputFile = sys.argv[1]
    else:
        sys.exit("Please give a hmmscan output file as the first argument")

    if len(sys.argv) > 3:
        eval = float(sys.argv[2])
        coverage = float(sys.argv[3])

    hmmscan_parse(inputFile, eval = eval, coverage = coverage)

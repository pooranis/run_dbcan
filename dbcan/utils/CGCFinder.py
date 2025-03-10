#!/usr/bin/env python3
##################################################
# CGCFinder v4
#
# Rewriten for dbCAN2
#
# Written by Le Huang under the supervision of
# Dr. Yin at NIU
#
#
# Last updated 12/24/18
#-updating info
#-adding stp
##################################################






# #open output file
# out = open(args.output, 'w+')

#global vars
cluster = [0, 0, 0, 0] #cazyme, tp, tf, stp
num_clusters = 0


#define boolean function to determine if a cluster meets cluster requirements
def isCluster(siggenes):
	global cluster
	if siggenes == 'all':
		if cluster[0] > 0 and cluster[1] > 0 and cluster[2] > 0 and cluster[3]:
			return True
	elif siggenes == 'tf':
		if cluster[0] > 0 and cluster[2] > 0:
			return True
	elif siggenes == 'tp':
		if cluster[0] > 0 and cluster[1] > 0:
			return True
	elif siggenes == 'stp':
		if cluster[0] > 0 and cluster[3] > 0:
			return True
	elif siggenes == 'tp+tf':
		if cluster[0] > 0 and cluster[1] > 0 and cluster[2]>0:
			return True
	elif siggenes == 'tp+stp':
		if cluster[0] > 0 and cluster[1] > 0 and cluster[3]>0:
			return True
	elif siggenes == 'tf+stp':
		if cluster[0] > 0 and cluster[2] > 0 and cluster[3]>0:
			return True
	else:
		print('Warning: invalid siggenes argument')
	return False

#define boolean function to detemine if a gene is important (a signature gene)
def isImportant(gene, siggenes):
	if gene == 'CAZyme':
		return True
	else:
		if gene == 'TC' and (siggenes == 'tp' or siggenes == 'all' or siggenes == 'tp+tf' or siggenes == 'tp+stp'):
			return True
		if gene == 'TF' and (siggenes == 'tf' or siggenes == 'all' or siggenes == 'tp+tf' or siggenes == 'tf+stp'):
			return True
		if gene == 'STP' and (siggenes == 'stp' or siggenes == 'all'or siggenes == 'tp+stp' or siggenes == 'tf+stp'):
			return True
	return False

def isSigGene(gene):
	if gene == 'CAZyme' or gene == 'TC' or gene == 'TF' or gene == 'STP':
		return True
	else:
		return False

#define function to increase the cluster count
def increaseClusterCount(gene):
	global cluster
	if gene == 'CAZyme':
		cluster[0] += 1
	elif gene == 'TC':
		cluster[1] += 1
	elif gene == 'TF':
		cluster[2] += 1
	elif gene == 'STP':
		cluster[3] += 1
	else:
		print("Warning: increaseClusterCount was called on bad functional domain")

#define function to search for a cluster once an important gene has been found
#this function also handles output
def startSearch(startRow, contig, distance, siggene, out):
	global cluster
	global num_clusters
	dis = distance
	index = startRow
	between = 0
	lastImportant = 0
	while index < len(contig):
		index += 1
		fd = contig[index][2]
		if isImportant(fd, siggene):
			increaseClusterCount(fd)
			lastImportant = index
			between = 0
		else:
			between += 1
		if between > dis or index >= (len(contig)-1):
			if isCluster(siggene):
				num_clusters += 1
				#output file columns
				#geneNumber type[2] downDis upDis CGC# contig[0] geneStart[3] geneEnd[4] geneID[8,ID] direc[6] note[8]
				for j in range(startRow, lastImportant + 1):
					fd = contig[j][2]
					if isSigGene(fd):
						upDown = findNear(contig, j, siggene)
						notes = contig[j][8].split(";")
						ID= ""
						for note in notes:
							if "ID" in note:
								ID = note.split("=")[1]
						row = [str(j), fd, str(upDown[1]), str(upDown[0]), 'CGC'+str(num_clusters), contig[j][0], contig[j][3], contig[j][4], ID, contig[j][6], contig[j][8]]
					else:
						row = [str(j), 'null', 'null', 'null', 'CGC'+str(num_clusters), contig[j][0], contig[j][3], contig[j][4], ID, contig[j][6]]
					try:
						row.append(contig[j][8])
					except:
						pass
					out.write('\t'.join(row) + '\n')
				out.write('+++++' + '\n')
			cluster = [0, 0, 0, 0]
			return index

#define function to find how close important genes are to each other
def findNear(contig, index, siggene):
	vals = ['null', 'null']
	k = index - 1
	l = index + 1
	while k >= 0:
		if isImportant(contig[k][2], siggene):
			vals[0] = index - k - 1
			break
		else:
			k -= 1
	while l <= len(contig) - 1:
		if isImportant(contig[l][2], siggene):
			vals[1] = l - index - 1
			break
		else:
			l += 1
	return vals


def cgc_finder(gffFile, distance, siggenes, output):
	global cluster
	global num_clusters

	#open output file
	out = open(output, 'w+')

	#load contig into an array
	contigs = {}
	with open(gffFile) as f:
		for line in f:
			row = line.rstrip().split('\t')
			if row[0] not in contigs:
				contigs[row[0]] = []
			contigs[row[0]].append(row)

	#loop through contig
	for key in contigs:
		contig = contigs[key]
		num_clusters = 0
		i = 0
		while i < len(contig) - 1:
			fd = contig[i][2]

			if isImportant(fd, siggenes):
				increaseClusterCount(fd)
				i = startSearch(i, contig, distance, siggenes, out)
			else:
				i += 1

	if output != 'none':
		out.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='CAZyme Gene Cluster Finder')

    parser.add_argument('gffFile', help='GFF file containing genome information')
    parser.add_argument('--distance', '-d', type=int, choices=[0,1,2,3,4,5,6,7,8,9,10], default=2, help='The distance allowed between two signature genes')
    parser.add_argument('--siggenes', '-s', choices=['all', 'tp', 'tf','stp','tp+tf','tp+stp','tf+stp'], default='all', help='Signature genes types required. all=CAZymes,TC,TF; tp=CAZymes,TC; tf=CAZymes,TF')
    parser.add_argument('--output', '-o', help='Output file name')

    args = parser.parse_args()
    cgc_finder(args.gffFile, args.distance, args.siggenes, args.output)

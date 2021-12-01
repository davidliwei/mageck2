import csv
import numpy as np

class Gene(object):
	def __init__(self,genename,chrno,start,end):
		self.genename = genename
		self.chrno = chrno
		self.start = int(start)
		self.end = int(end)

class sgRNA(object):
	def __init__(self,sgrna,genename,chrno,start,end):
		self.sgrna = sgrna
		self.genename = genename
		self.chrno = chrno
		self.start = int(start)
		self.end = int(end)	
		
def readGenePosBED(pos_file):
	'''takes as input a tab-delimited file in BED format containing data on the 
	chromosomal position of genes and returns a dictionary of dictionary of objects 
	containing chromosome number, start & end positions for each gene indexed by gene 
	and indexed by chromosome number for genes represented in the inputted file'''

	with open(pos_file,'r') as f:
		reader = csv.reader(f,delimiter='\t')
		header = reader.next()

		geneIdx = 3
		chrnoIdx = 0
		startIdx = 1
		endIdx = 2

		pos_dict = {}
		for line in reader:
			if line[chrnoIdx] not in pos_dict:
				pos_dict[line[chrnoIdx]] = {}
			pos_dict[line[chrnoIdx]][line[geneIdx]] = \
				Gene(line[geneIdx],line[chrnoIdx],line[startIdx],line[endIdx])

	return pos_dict


def gene2sgrnaPos(genepos_dict,sgrna2genelist):
	'''takes in a dictionary of Gene objects indexed by gene symbol/name and returns
	a corresponding dictionary of dictionaries of sgRNA objects indexed by sgRNA name
	indexed by chromosome number (as was done with the Gene objects dictionary)'''

	# invert sgrna2genelist to get all sgRNAs associated with each gene
	gene2sgrna_dict = {gene: [] for gene in sgrna2genelist.values()}
	[gene2sgrna_dict[sgrna2genelist[sgrna]].append(sgrna) for sgrna in sgrna2genelist]

	sgrnapos_dict = {chrno: {} for chrno in genepos_dict}
	for chrno in genepos_dict:
		for gene in genepos_dict[chrno]:
			if gene in gene2sgrna_dict:
				for sgrna in gene2sgrna_dict[gene]:
					start = genepos_dict[chrno][gene].start
					end = genepos_dict[chrno][gene].end
					sgrnapos_dict[chrno][sgrna] = sgRNA(sgrna,gene,chrno,start,end)

	return sgrnapos_dict
	
def getGeneWindowScores(pos_dict,chrno,score_dict,window_size,step_size,summary_op=np.mean):
	'''segments a chromosome into windows of an inputted window size, performs
	a summary operation (i.e. mean) on the scores corresponding to genes in each
	window, and returns a list of these scores; also returns a list of gene lists
	with indices corresponding to the outputted window scores'''

	# sort gene objects from the inputted chromosome number by start position
	sortedGenes_list = sorted(pos_dict[chrno].values(),key=lambda x: x.start)

	# initialize list of gene lists
	genesWindow_list = []

	# get last position in chromosome
	lastPos = sortedGenes_list[len(sortedGenes_list)-1].end

	firstGeneIdx = 0
	windowSummary_list = []
	currentEnd = 0
	for window_end in range(window_size,lastPos,step_size):

		window_scores = [] # stores
		genesInWindow = [] # stores genes in the current window
		GeneIdx = firstGeneIdx
		while currentEnd < window_end:
			# mark index of first gene to be found in window
			if len(window_scores) == 1:
				firstGeneIdx = GeneIdx

			currentGenename = sortedGenes_list[GeneIdx].sgrna # .genename
			currentStart = sortedGenes_list[GeneIdx].start
			currentEnd = sortedGenes_list[GeneIdx].end
			GeneIdx += 1

			# get score associated with gene in window
			if currentStart > window_end-window_size and currentGenename in score_dict:
				window_scores.append(score_dict[currentGenename])
				genesInWindow.append(currentGenename)

		# get summary statistic for scores of genes within window
		# (if there are more than 4 genes/scores within a window)
		window_scores = np.array(window_scores)
		if len(window_scores) > 4:
			windowSummary_list.append(summary_op(window_scores))
			genesWindow_list.append(genesInWindow)

	return windowSummary_list, genesWindow_list

def allWindowScoresGenes(pos_dict,score_dict,window_size,step_size,summary_op=np.mean):
	'''computes summary scores for all windows in each chromosome and outputs
	a list of all of the summary scores along with a corresponding list of gene
	lists associated with the windows'''

	allWindowScores = []
	allWindowGenes = []
	for chrno in pos_dict.keys():
		summary_scores, genesWindow_list = getGeneWindowScores(pos_dict,chrno, \
			score_dict,window_size,step_size,summary_op)
		allWindowScores.extend(summary_scores)
		allWindowGenes.extend(genesWindow_list)

	return allWindowScores, allWindowGenes

def geneWindowAvgScore_dict(allWindowGenes,allWindowScores):
	'''takes in a list of gene lists (each of which corresponds to the set of
	genes positioned within a window) and returns a dictionary of averaged 
	window scores for all windows associated with each gene (i.e. get average 
	window scores for each gene)'''

	# get indices of windows that contain each gene
	gene_dict = {}
	for i in range(len(allWindowGenes)):
		for gene in allWindowGenes[i]:
			if gene not in gene_dict:
				gene_dict[gene] = []
			gene_dict[gene].append(i)

	# determine average window scores for each gene
	return {gene: np.mean([allWindowScores[scoreIdx] for scoreIdx in \
		gene_dict[gene]]) for gene in gene_dict}	

def calcLFC(tabctrl,tabtest):
	'''calculates log-fold change associated with each sgRNA and returns
	a dictionary of log-fold change values indexed by sgRNA name'''

	return {sgrna: np.log2(np.mean(tabtest[sgrna])/np.mean(tabctrl[sgrna])) \
		for sgrna in tabctrl if np.mean(tabctrl[sgrna]) > 0 \
		and np.mean(tabtest[sgrna]) > 0}

def mapCNVtoGene(CNV_dict,sgrna2genelist):
	'''takes a dictionary of CNV estimates indexed by a Cytoband object and
	returns a dictionary of CNV estimates for each gene indexed by gene'''

	geneCNV_dict = {}
	for region in CNV_dict:
		for sgrna in region.sgrnas:
			if sgrna2genelist[sgrna] not in geneCNV_dict:
				geneCNV_dict[sgrna2genelist[sgrna]] = CNV_dict[region]

	return geneCNV_dict

def writeCNVestimates2file(geneCNV_dict_list,output_prefix,cell_line_list):
	'''writes the CNV estimates from multiple experiments (analyzed using MAGeCK-MLE)
	to a tab-delimited .txt file'''

	f = open(output_prefix + '.CNVestimates.txt','wb')
	writer = csv.writer(f,delimiter='\t')
	writer.writerow(['SYMBOL']+cell_line_list)
	genelist = sorted(geneCNV_dict_list[0].keys())
	[writer.writerow([gene] + [str(geneCNV_dict[gene]) for geneCNV_dict in \
		geneCNV_dict_list]) for gene in genelist]
	f.close()

def estimateCNV(genescore_dict,allWindowScores):
	'''estimates CNV for each gene from the scores associated with the
	corresponding genes'''

	geneCNV_dict = {}
	scoresSD = np.nanstd(allWindowScores)
	CNVestMean = -1*np.nanmean(allWindowScores)

	for gene in genescore_dict:
		score = genescore_dict[gene]
		CNVest = 1./scoresSD*(-1*score-CNVestMean)+CNVestMean
		geneCNV_dict[gene] = np.clip(CNVest,-4,4)

	return geneCNV_dict

def mageckCNVestimation(genepos_file,tabctrl,tabtest,sgrna2genelist,cell_line,output_prefix,\
	window_size=2000000,step_size=100000,filterThreshold=1.):
	'''unified function for estimating CNV given a file with gene position
	information and control/test read counts (for MAGeCK)'''

	# read gene position data and map to sgRNA positions
	genepos_dict = readGenePosBED(genepos_file)
	sgrnapos_dict = gene2sgrnaPos(genepos_dict,sgrna2genelist)

	# calculate LFC values from read count tables
	score_dict = calcLFC(tabctrl,tabtest)

	# calculate window LFC scores and determine sgRNAs in each window
	allWindowScores, sgrnasWindow = allWindowScoresGenes(sgrnapos_dict,score_dict,\
		window_size,step_size,summary_op=np.mean)
	sgrnaWindowScores = geneWindowAvgScore_dict(sgrnasWindow,allWindowScores)

	# get scores associated with each gene
	geneLFC_dict = {sgrna2genelist[sgrna]: score for (sgrna,score) in sgrnaWindowScores.items()}

	# estimate CNV from LFC scores
	geneCNV_dict = estimateCNV(geneLFC_dict,allWindowScores)

	# write CNV estimates to file
	writeCNVestimates2file([geneCNV_dict],output_prefix,[cell_line])

##########

def ctrl_test_from_designmatrix(cttab,designmat):
	'''returns a list of control and test read count tables (stored in tuples of 
	dictionaries) corresponding to the control/test pairings and experiments 
	specified in the design matrix'''

	ctrltest_dicts = []

	# get indices of design matrix corresponding to control samples
	designRowSums = np.sum(designmat,axis=1)
	ctrlInds = np.array([i for i in range(designmat.shape[0]) if designRowSums[i]==1])

	ctrltab = {key: list(np.array(val)[ctrlInds]) for (key,val) in cttab.items()}

	for i in range(designmat.shape[1]-1):
		testInds = [j for j in range(designmat.shape[0]) if designmat[j,i+1]==1]
		testtab = {key: list(np.array(val)[testInds]) for (key,val) in cttab.items()}
		ctrltest_dicts.append((ctrltab,testtab))

	return ctrltest_dicts

def mageckmleCNVestimation(genepos_datafile,cttab,designmat,sgrna2genelist,beta_labels,output_prefix,\
	window_size=2000000,step_size=100000):
	'''unified function for estimating CNV given a file with gene position
	information and control/test read counts (for MAGeCK-MLE)'''

	# read gene position data and map to sgRNA positions
	genepos_dict = readGenePosBED(genepos_datafile)
	sgrnapos_dict = gene2sgrnaPos(genepos_dict,sgrna2genelist)

	# separate control and test read count tables
	ctrltest_dicts = ctrl_test_from_designmatrix(cttab,designmat)

	geneCNVestimates = []
	for (tabctrl,tabtest) in ctrltest_dicts:

		# calculate LFC values from read count tables
		score_dict = calcLFC(tabctrl,tabtest)

		# calculate window LFC scores and determine sgRNAs in each window
		allWindowScores, sgrnasWindow = allWindowScoresGenes(sgrnapos_dict,score_dict,\
			window_size,step_size,summary_op=np.mean)
		sgrnaWindowScores = geneWindowAvgScore_dict(sgrnasWindow,allWindowScores)

		# get scores associated with each gene
		geneLFC_dict = {sgrna2genelist[sgrna]: score for (sgrna,score) in sgrnaWindowScores.items()}

		# estimate CNV from LFC scores
		geneCNV_dict = estimateCNV(geneLFC_dict,allWindowScores)

		geneCNVestimates.append(geneCNV_dict)
	
	# write CNV estimates to file
	writeCNVestimates2file(geneCNVestimates,output_prefix,beta_labels) 


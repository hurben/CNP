def testset():

	with open('//data/project/hurben/dataset/Per2/PER2_WAT_KO.DEG.list') as f:
	    lines = f.read().splitlines()

	KO_gene = 'Per2'
	context = ['Lipid Metabolism','inflammatory response', 'adipocyte differentiation', 'adipogenesis']

	return lines, KO_gene, context


def Text_to_List(input_file):

	DEG_list = []

	open_input_file = open(input_file,'r')
	input_file_readlines = open_input_file.readlines()
	
	for i in range(len(input_file_readlines)):
		read = input_file_readlines[i]
		read = read.replace('\n','')
		
		#check whether there are any duplicate genes
		if read in DEG_list:
			print ("[WARNING] Duplicate genes in Input file's entry. Caution is required")
			print ("[WARNING] Checking the whole process is adviced")

		if read not in DEG_list:
			DEG_list.append(read)

	return DEG_list


def attain_seed_gene_list(user_context, KO_gene, DEG_list):

	print ('[Notice] Proceeding BEST')
	
	bestQuery = best.BESTQuery({"keywordA":[user_context], "keywordB":[KO_gene], "filterObjectName":"","topN":100})
	relevantEntities = best.getRelevantBioEntities(bestQuery)
	print ((relevantEntities['genes']))
	print (len(relevantEntities['genes']))

	LEN_seed_genes = len(relevantEntities['genes'])
	seed_gene_list = []

	for i in range(LEN_seed_genes):
		seed_gene = relevantEntities['genes'][i]['entityName']
		seed_gene = ("%s%s" % (seed_gene[0].upper(), seed_gene[1:].lower()))

		if seed_gene not in seed_gene_list:
			if seed_gene not in DEG_list:
				seed_gene_list.append(seed_gene)
		else:
			print ("some error #1")

	return seed_gene_list


def StringDB_to_dict():
	print ('[FL_MLV] START StringDB_to_dict')

	StringDB_dir = '//data/project/hurben/database/MED_STRING_PPI_Mus_musculus_Symbol.txt'
#	StringDB_dir = '//data/project/hurben/database/HIGH_STRING_PPI_Mus_musculus_Symbol.txt'
	StringDB_open = open(StringDB_dir,'r')
	StringDB_readlines = StringDB_open.readlines()

	StringDB_dict = {}
	Gene_list = []

	for i in range(len(StringDB_readlines)):
		read = StringDB_readlines[i]
		read = read.replace('\n','')
		token = read.split('\t')
		node_1 = token[0]
		node_2 = token[1]
		#if node_1 not in StringDB_dict.keys():
		try :
			StringDB_dict[node_1].append(node_2)
		except KeyError:
			StringDB_dict[node_1] = [node_2]

		Gene_list.append(node_1)
		Gene_list.append(node_2)

	Gene_list = list(set(Gene_list))
	SortGene_list = sorted(Gene_list)
	print ("[FL_MLV] Total # of Genes: ", len(SortGene_list))
	print ("[FL_MLV] Total # of StringDB_dict Keys: ", len(StringDB_dict.keys()))
	print ('[FL_MLV] END StringDB_to_dict')

	return StringDB_dict, SortGene_list



def Create_RWR_Condition_AdjacencyMatrix_part1(StringDB_dict, total_gene_list):
#Express condition should be unified.
#ex) only UP DEG or only DOWN DEG
	print ('[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : START')
	topology_dict = {}

	for node1 in StringDB_dict.keys():
		if node1 in total_gene_list:
		#[1] If PPI node 1 is DEG gene

			topology_dict[node1] = []
			node2_list = StringDB_dict[node1]

			for node2 in node2_list:
				if node2 in total_gene_list:
					#[2] And PPI node 2 is DEG gene, append
					topology_dict[node1].append(node2)

	print ("[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : END")

	##debug
	debug_txt = open('debug.txt','w')
	for key in topology_dict.keys():
		debug_txt.write('%s\n' % key)
		dict_gene_list = topology_dict[key]
		debug_txt.write('\n'.join(dict_gene_list))


	return topology_dict
	#DEGtopology_dict : contains every node-edges of every conditions.
	#total_gene_list : every genes of every conditions


def Create_RWR_Condition_AdjacencyMatrix_part2(topology_dict, total_gene_list, tag):
	print ('[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 2 : START')
	adj_matrix_file_name = str(tag) + '.adj.matrix'
	adj_matrix_txt = open(adj_matrix_file_name,'w')

	print ('[FL_stringDB_RWR] Number of genes for Adjacency matrix : ', len(total_gene_list))
	for node1 in total_gene_list:
		if node1 in topology_dict.keys():
			adj_matrix_txt.write('\t' + str(node1))
	adj_matrix_txt.write('\n')

	for node1 in total_gene_list:
		if node1 in topology_dict.keys():

			adj_matrix_txt.write(str(node1))
			#[1] write row name

			for node2 in total_gene_list:
				if node2 in topology_dict.keys():
					if node1 == node2:
					#[2-1] check current position. is duped
						adj_matrix_txt.write('\t0')

					if node1 != node2:
					#[2-2] if not duped,
						node2_list = topology_dict[node1]

						if node2 in node2_list:
							adj_matrix_txt.write('\t1')
						if node2 not in node2_list:
							adj_matrix_txt.write('\t0')
			adj_matrix_txt.write('\n')
	adj_matrix_txt.close()


def Create_RWR_p0_vector(topology_dict, total_gene_list, seed_gene_list, tag): 

	p0_vector_txt = open(str(tag) + '.p0.vector', 'w')

	for gene in total_gene_list:
		if gene in topology_dict.keys():
			p0_vector_txt.write(str(gene))

			#[2] if gene is in intersection, it is considered as seed = p0
			if gene in seed_gene_list:
				p0_vector_txt.write('\t1\n')

			#[2] if gene is not in intersection, it is not considered as seed = p0
			if gene not in seed_gene_list:
				p0_vector_txt.write('\t0\n')



def Run_RWR(tag):

	path_to_rwr_script = '//data/project/hurben/program/RWR.R'
	p0_vector_file = str(tag) + '.p0.vector'
	adj_matrix_file = str(tag) + '.adj.matrix'
	rwr_result = str(tag) + '.rwr'
	Rscript_RWR_cmd = 'Rscript ' + str(path_to_rwr_script) + ' ' + str(adj_matrix_file) + ' ' + str(p0_vector_file) + ' ' + str(rwr_result)
	print (Rscript_RWR_cmd)
	os.system(Rscript_RWR_cmd)


def organize_rwr_results(tag, rwr_summary_dict, available_deg_list):
	print (tag)

	rwr_file = str(tag) + '.rwr'
	rwr_open = open(rwr_file,'r')

	rwr_readlines = rwr_open.readlines()
	score_dict = {}

	for i in range(1, len(rwr_readlines)):
		read = rwr_readlines[i]
		read = read.replace('\n','')
		token = read.split()

		try :
			gene = token[0]
			prob = float(token[1])
		except IndexError:
			print ('############################ERROR')
			print (read)

		score_dict[gene] = prob


	sorted_score_list = sorted(score_dict.items(), key=operator.itemgetter(1), reverse=True)
	for i in range(len(sorted_score_list)):
		gene = sorted_score_list[i][0]
		score = sorted_score_list[i][1]

		if gene in DEG_list:
			rwr_summary_dict[gene,tag] = score
			if gene not in available_deg_list:
				available_deg_list.append(gene)

	return rwr_summary_dict, available_deg_list


def rwr_summary(rwr_summary_dict, available_deg_list, context_list):

#	import numpy as np
#	for gene in DEG_list:
#		score_list = []
#		for i in range(20):
#			score = float(rwr_summary_dict[gene,i])
#			score_list.append(score)
#		std = np.std(score_list)
#		mean = np.mean(score_list)
	print ('#------ Summarizing ------- ')
	print ('Total considered candidates : %s' % len(available_deg_list))


	f = open('rwr.summary.table.txt','w')
	for context in context_list:
		f.write('\t' + str(context))

	for i in range(iteration):
		f.write('\t' + str(i))
	f.write('\n')

	for gene in available_deg_list:
		f.write(gene)

		for context in context_list:
			context = context.replace(' ','_')
			score = rwr_summary_dict[gene, context]
			f.write('\t' + str(score))

		for i in range(iteration):
			score = rwr_summary_dict[gene,i]
			f.write('\t' + str(score))

		f.write('\n')


def random_seed_gene(ppi_dict, DEG_list, average_seed_gene):

	seed_gene_list = []
	LEN_ppi_dict = len(ppi_dict.keys())

	for i in range(round(average_seed_gene)):

		gene = list(ppi_dict.keys())[randrange(LEN_ppi_dict)]

		if gene not in seed_gene_list:
			if gene not in DEG_list:
				seed_gene_list.append(gene)

	return seed_gene_list
		


if __name__ == '__main__':

	import sys
	import argparse
	import subprocess
	import os
	import operator
	from random import randrange
	
	#importing Korea univ. boss system
	program_dir = os.path.split(os.path.realpath('__file__'))
	api_lib_path = os.path.abspath('//data/project/hurben/BEST/dmis')
	sys.path.insert(0, str(api_lib_path))
	import best

	#initialize
	GO_term_list = Text_to_List('//data/project/hurben/database/GO_term.list')
	ppi_dict, ppi_gene_list = StringDB_to_dict()
	rwr_summary_dict = {}
	available_deg_list = []

		
	DEG_list, KO_gene, context_list = testset()
	iteration = 100
	seed_gene_count = 0

	print ("Number of DEG genes", len(DEG_list))
#	print ("Number of seed genes", len(seed_gene_list))

	for context in context_list:

		print (context)
		seed_gene_list = attain_seed_gene_list(context, KO_gene, DEG_list)
		seed_gene_count += len(seed_gene_list)
		total_gene_list = DEG_list + seed_gene_list
		tag = context
		tag = tag.replace(' ','_')

		topology_dict = Create_RWR_Condition_AdjacencyMatrix_part1(ppi_dict, total_gene_list)
		Create_RWR_Condition_AdjacencyMatrix_part2(topology_dict, total_gene_list, tag)
		Create_RWR_p0_vector(topology_dict, total_gene_list, seed_gene_list, tag)
		Run_RWR(tag)
		rwr_summary_dict, available_deg_list = organize_rwr_results(tag,rwr_summary_dict,available_deg_list)
	
	average_seed_gene = seed_gene_count / 4
	print (average_seed_gene)

	for i in range(iteration):

		seed_gene_list = random_seed_gene(ppi_dict, DEG_list, average_seed_gene)
		total_gene_list = DEG_list + seed_gene_list
		tag = i

		topology_dict = Create_RWR_Condition_AdjacencyMatrix_part1(ppi_dict, total_gene_list)
		Create_RWR_Condition_AdjacencyMatrix_part2(topology_dict, total_gene_list, tag)
		Create_RWR_p0_vector(topology_dict, total_gene_list, seed_gene_list, tag)
		Run_RWR(tag)
		rwr_summary_dict, available_deg_list = organize_rwr_results(tag, rwr_summary_dict,available_deg_list)

	rwr_summary(rwr_summary_dict, available_deg_list, context_list)



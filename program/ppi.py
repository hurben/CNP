#PROGRAM DESIGNED TO BE RUNNED by "Analysis_Initialize"
#VER DEG
import sys
import os
import os.path
import argparse
import FL_MLV_Methods



def REMOVE_FC(gene):
	gene = gene.split('_')[0]
	return gene

def REMOVE_FC_FROM_LIST(gene_fc_list):
	
	gene_list = []
	for gene_fc in gene_fc_list:
		gene = REMOVE_FC(gene_fc)
		gene_list.append(gene)
	return gene_list

def INPUT_CHECK(gene_list_profile, method, output):

	if gene_list_profile == None:
		print "Gene profile not given. Terminating Program"
		quit()

	if method == None:
		print "Method not given. Terminating Program"
		quit()

	if output == None:
		output = 'default'

	return output

def CHECK_UNABLE_GENES(condition_dict, ppi_gene_list, result_file_string):

	#CHECK unable genes
	unable_gene_list = []
	for condition in condition_dict.keys():
		deg_genes = condition_dict[condition]

		for gene in deg_genes:
			gene = REMOVE_FC(gene) #NEED
			if gene not in ppi_gene_list:
				unable_gene_list.append(gene)

	unable_gene_list = list(set(unable_gene_list))

	open_error_txt = open(str(result_file_string) + '.unable.degs.txt','w')
	for gene in unable_gene_list:
		open_error_txt.write(str(gene) + '\n')
	open_error_txt.close()

	return unable_gene_list


def MAIN(condition_dict, method, cut_off, cut_off_rank,selected_condition, interest_condition, gene_condition_dict):

	result_dict = {}

	seed_gene_list, seed_condition_list = FL_MLV_Methods.Venndianet_RWR().Set_Seed_Genes(condition_dict, selected_condition)
	print 'Seeds are given :', selected_condition, len(seed_gene_list)

	#FAST NETWORK VERSION
	if method == 'fast':
		total_gene_list = FL_MLV_Methods.Venndianet_RWR().Set_subnetwork_Genes(condition_dict, selected_condition, interest_condition)

		DEG_topology_dict, common_gene_list = FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_subnetwork_part1(total_gene_list, ppi_dict, unable_gene_list)
		FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_subnetwork_part2(DEG_topology_dict, common_gene_list, result_file_string, unable_gene_list)

	#TOTAL NETWORK VERSION
	if method != 'fast':
		DEG_topology_dict, common_gene_list = FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_part1(condition_dict, ppi_dict, unable_gene_list)
		FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_part2(DEG_topology_dict, common_gene_list, result_file_string, unable_gene_list)


	FL_MLV_Methods.Venndianet_RWR().Create_RWR_p0_vector(DEG_topology_dict, seed_gene_list, result_file_string)
	FL_MLV_Methods.Venndianet_RWR().Run_RWR(result_file_string)

	if cut_off != None:
		rwr_cut_dict, gene_prob_score_dict = FL_MLV_Methods.Venndianet_RWR().Cut_RWR(result_file_string, gene_condition_dict, cut_off, seed_condition_list)
	if cut_off_rank != None:
		rwr_cut_dict, gene_prob_score_dict = FL_MLV_Methods.Venndianet_RWR().Cut_RWR_by_Rank(result_file_string, gene_condition_dict, cut_off_rank, seed_condition_list)

	return rwr_cut_dict, gene_prob_score_dict


def CREATE_OUTPUT(rwr_cut_dict, gene_prob_score_dict, result_file_string, condition_of_interest,cut_off, cut_off_rank):

	if cut_off == None:
		cut_off = cut_off_rank

	condition_of_interest_list = list(map(int, condition_of_interest.split('+')))

	output_txt = open(str(result_file_string) + '.rwr.' + str(cut_off) + '.cut','w')

	##### 1 ######
	output_table_txt = open(str(result_file_string) + '.rwr.' + str(cut_off) + '.cut.table.forRplot','w')
	output_table_txt.write('condition\tprob\trank\n')

	conditionID_list = rwr_cut_dict.keys()
	
	import operator
	import scipy.stats as SS

	prob_rank_dict = {}
	prob_list = []

	sorted_gene_prob_score_list = sorted(gene_prob_score_dict.items(), key=operator.itemgetter(1), reverse=True)

	for sorted_item in sorted_gene_prob_score_list:
		gene = sorted_item[0]

		prob = float(gene_prob_score_dict[gene])
		prob_list.append(prob)

	prob_list = list(set(prob_list))
	ranked_prob_list = len(prob_list) - SS.rankdata(prob_list)

	for i in range(len(prob_list)):
		prob = float(prob_list[i])
		rank = float(ranked_prob_list[i])
		prob_rank_dict[prob] = rank

	for sorted_item in sorted_gene_prob_score_list:
		gene = sorted_item[0]

		for conditionID in conditionID_list:

			if conditionID != "seed":
				if conditionID in condition_of_interest_list:
					prob = gene_prob_score_dict[gene]
					rank = prob_rank_dict[prob]

					if gene in rwr_cut_dict[conditionID]:
						output_table_txt.write('C' + str(conditionID) + '\t' + str(prob) + '\t' + str(rank) + '\n')
	output_table_txt.close()


	##### 2 ######

	output_table_txt = open(str(result_file_string) + '.rwr.' + str(cut_off) + '.cut.forRstats','w')
	output_table_txt.write('gene')
	for conditionID in rwr_cut_dict.keys():
		output_table_txt.write('\t' + str(conditionID))
	output_table_txt.write('\n')

	conditionID_list = rwr_cut_dict.keys()

	for gene in gene_prob_score_dict.keys():
		output_table_txt.write(str(gene))

		for conditionID in conditionID_list:

			prob = gene_prob_score_dict[gene]
			##
			prob = prob_rank_dict[prob]
			##
			if gene in rwr_cut_dict[conditionID]:
				output_table_txt.write('\t' + str(prob))
			if gene not in rwr_cut_dict[conditionID]:
				output_table_txt.write('\tNull')
		output_table_txt.write('\n')

	output_table_txt.close()

	output_txt.write('gene\tprob\n')
	for conditionID in conditionID_list:

		for gene in rwr_cut_dict[conditionID]:
			prob = gene_prob_score_dict[gene]
			output_txt.write(str(gene) + '\t' + str(prob) + '\n')

	output_txt.close()
		


if __name__ == '__main__':

	FL_MLV_Methods.WELCOME_MESSAGE()
	print '[VENN_PPI] START'
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', dest = 'gene_list')

	parser.add_argument('-c', '--cutoff', dest = 'cutoff')
	parser.add_argument('-cr', '--cutoff_rank', dest = 'cutoff_rank')

	parser.add_argument('-o', '--output', dest = 'result_file')
	parser.add_argument('-m', '--method', dest = 'method')

	parser.add_argument('-s', '--seed', dest = 'seed_condition')
	parser.add_argument('-ci', '--condition_of_interest', dest = 'condition_of_interest')
	parser.add_argument('-ppi', '--ppi_method', dest = 'ppi_method')

	args = parser.parse_args()
	gene_list_profile = args.gene_list
	cut_off = args.cutoff
	cut_off_rank = args.cutoff_rank
	ppi = args.ppi_method
	result_file_string = args.result_file
	method = args.method
	seed_condition = args.seed_condition
	interest_condition = args.condition_of_interest


	if cut_off != None:
		print "CUTOFF: ", cut_off
	if cut_off_rank != None:
		print "CUTOFF RANK: ", cut_off_rank

	if method == None:
		method = 'slow'
	
	
	gene_list_file = open(gene_list_profile,'r')

	if ppi == None or ppi == 'med':
		ppi_profile_file = open('//var/www/htdocs/MLV/analysis/ppi/MED_STRING_PPI_Mus_musculus_Symbol.txt','r')
	if ppi == 'high':
		ppi_profile_file = open('//var/www/htdocs/MLV/analysis/ppi/HIGH_STRING_PPI_Mus_musculus_Symbol.txt','r')

	gene_list_file_readlines = gene_list_file.readlines()
	ppi_profile_file_readlines = ppi_profile_file.readlines()

	#FL----
	condition_dict, gene_list, gene_condition_dict = FL_MLV_Methods.MLV_Profile_to_dict(gene_list_file_readlines)
	ppi_dict, ppi_gene_list = FL_MLV_Methods.StringDB_to_dict(ppi)
	#FL----

	unable_gene_list = CHECK_UNABLE_GENES(condition_dict, ppi_gene_list, result_file_string)

	rwr_cut_dict, gene_prob_score_dict = MAIN(condition_dict, method, cut_off, cut_off_rank, seed_condition, interest_condition, gene_condition_dict)
	CREATE_OUTPUT(rwr_cut_dict, gene_prob_score_dict,result_file_string, interest_condition, cut_off, cut_off_rank)

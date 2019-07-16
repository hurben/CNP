################################################################
#api.py													15.02.09 -> 15.07.28 -> 19.07.16  Initiating PLAN:B
#scripted in python 3.4.2
#
#Description : ver max
#
#
#cmd : python api.py -i [input]  -c [context] -o [output]
################################################################
#input
#string1\n
#string2\n
#string3\n

#context
#string

#output
#string


#CHECKING OPTIONS
#checking whether input are correctly given
def Check_Options(input_file, context_aware, output_file):
	if input_file != None:
		print ('[Notice] Reading Input file')
		caff_gene_results = Text_to_List(input_file)

	if input_file == None:
		print ('[Error] Input file not given')
		print ('[Error] terminating program')
		quit()

	if context_aware == None:
		print ('[Notice] Context : No special interest')
	if context_aware != None:
		print ('[Notice] Context : ' + str(context_aware))
		context_aware = [str(context_aware)]
		
	if output_file == None:
		output_file = 'Triple_filter.final.result'

	return caff_gene_results, context_aware, output_file




#Without context consideraton	
def Boss_to_Dict(caff_gene_results, ko_gene):

	result_dict = {}
	print ('[Notice] Proceeding BEST')
	
		
	for gene in caff_gene_filter_results:
		query = [gene]
		bossQuery = boss.BOSSQuery({"keywordA":[ko_gene], "keywordB":query, "filterObjectName":"","topN":1})
		relevantEntities = boss.getRelevantBioEntities(bossQuery)


		try : 
			result_dict[gene] = relevantEntities['genes'][0]['score']
		except KeyError: 
			print ("[Notice] BEST cannot find " +str(gene) + ". It will be excluded from candidates.")
		#print (relevantEntities['genes'][0]['score']) #indicates rank. and dict inside dict

	return result_dict
	

#With context consideraton	
def Boss_to_Dict_with_Context(context, caff_gene_filter_results):

	result_dict = {}
	print ('[Notice] Proceeding BEST')
		
	for gene in caff_gene_filter_results:
		query = [gene]
		bossQuery = boss.BOSSQuery({"keywordA":context, "keywordB":query, "filterObjectName":"", "topN":1})
#		bossQuery = boss.BOSSQuery({"keywordA":["cancer"], "keywordB":["BRCA1", "EGFR"], "topN":1})
		relevantEntities = boss.getRelevantBioEntities(bossQuery)

		try : 
			result_dict[gene] = relevantEntities['genes'][0]['score']
		except KeyError: 
			print ("[Notice] BEST cannot find " +str(gene) + ". It will be excluded from candidates.")
		#print (relevantEntities['genes'][0]['score']) #indicates rank. and dict inside dict

	return result_dict
	




#Text file should have a format of 1 column
#ex) gene1\ngene2\ngene3\ngene4
def Text_to_List(input_file):

	list = []

	open_input_file = open(input_file,'r')
	input_file_readlines = open_input_file.readlines()
	
	for i in range(len(input_file_readlines)):
		read = input_file_readlines[i]
		read = read.replace('\n','')
		
		#check whether there are any duplicate genes
		if read in list:
			print ("[WARNING] Duplicate genes in Input file's entry. Caution is required")
			print ("[WARNING] Checking the whole process is adviced")

		if read not in list:
			list.append(read)

	return list

	
#Result making	
def	Result_Making(output_file, sorted_result_dict, result_dict):


	result_text = open(str(output_file),'w')

	result_text.write('RANK\tGENE\tSCORE\n')

	for i in range(len(sorted_result_dict)):
		gene = sorted_result_dict[i][0]
		print (gene)
		score = result_dict[gene]

		result_text.write(str(i) + '\t' + str(gene) + '\t' + str(score) +'\n')
		


if __name__ == '__main__':

	import sys
	import argparse
	import subprocess
	import os


	#importing Korea univ. boss system
	program_dir = os.path.split(os.path.realpath(__file__))
	api_lib_path = os.path.abspath(str(program_dir[0])+'/BEST/dmis')
	sys.path.insert(0, str(api_lib_path))
#	import boss
	import best

	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input', dest = 'input_file', help='List of genes resulted from "Triple_filter.py"')
	parser.add_argument('-c','--context', dest = 'context_aware', help='Context that you are interested"')
	parser.add_argument('-ko','--knockout', dest = 'knockout_gene', help='Knockout gene symbol"')
	parser.add_argument('-o','--output', dest = 'output_file', help='Define output file name')

	args= parser.parse_args()

	input_file = args.input_file
	output_file = args.output_file
	context_aware = args.context_aware
	ko_gene = args.knockout_gene


	print ('------------------------------------')
	print (': ::      CAFF-GENE & BEST     ::  :')
	print (':                                  :')
	print (':     >>   Result calling   <<     :')
	print (':                                  :')
	print (':                            1.0.0 :')
	print ('------------------------------------')


	if ko_gene == None:
		print ("Please give KO gene symbol")
		quit()

	if context_aware == None:
		context_aware = 'null'
		print ('null')

	if '_' in context_aware:

		context_aware = context_aware.replace('_',' ')
		context_aware = '"' + str(context_aware) + '"'




	caff_gene_filter_results, context_aware, output_file = Check_Options(input_file, context_aware, output_file)
	result_dict = Boss_to_Dict(caff_gene_filter_results, ko_gene)

	if context_aware[0] != 'null':

		result_context_dict = Boss_to_Dict_with_Context(context_aware, caff_gene_filter_results)

		gene_list = []

		for gene in result_context_dict.keys():
			gene_list.append(gene)

		for gene in result_dict.keys():
			gene_list.append(gene)

		print ('Total genes', len(gene_list))
		uniq_gene_list = list(set(gene_list))
		print ('uniq genes', len(uniq_gene_list))


		final_dict = {}

		for gene in uniq_gene_list:

			if gene in result_context_dict.keys() and gene in result_dict.keys():
				
				if float(result_context_dict[gene]) >= float(result_dict[gene]):
			
					final_dict[gene] = result_context_dict[gene]
					print ('both',gene, result_dict[gene],result_context_dict[gene])

				if float(result_context_dict[gene]) <= float(result_dict[gene]):

					final_dict[gene] = result_dict[gene]
					print ('both',gene, result_context_dict[gene])

			if gene in result_context_dict.keys() and gene not in result_dict.keys():

				final_dict[gene] = result_context_dict[gene]
				print ('context',gene, result_context_dict[gene])

			if gene not in result_context_dict.keys() and gene in result_dict.keys():

				final_dict[gene] = result_dict[gene]
				print ('ko',gene, result_dict[gene])



		print ('Total KEY size : ',len(final_dict.keys()))
	



		sorted_final_dict = sorted(final_dict.items(), key=lambda x: x[1], reverse=True)
		Result_Making(str(output_file) +'.context', sorted_final_dict, final_dict)




	else:
		sorted_result_dict = sorted(result_dict.items(), key=lambda x: x[1], reverse=True)

		Result_Making(str(output_file)+'.noncontext', sorted_result_dict, result_dict)


		
	


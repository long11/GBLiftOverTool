from parser_v5 import *

def initialize_file(filename, filename1_size, filename2_size, graph):
	"""
	Initialize the graph from files

	Args:
		filename 	   (string): Filename
		filename1_size (string): filename1 chromosome size
		filename2_size (string): filename2 chromosome size
		graph 	  	    (class): Class graph
	"""

	extracted, chrom1_size, chrom2_size = parse_file(filename, filename1_size, 
													 filename2_size, graph)

	#First chromosome
	unique_id, genome1_start_index_graph, length_common, length_unique_genome1 = \
			implement_first_chromosome_to_graph(extracted, chrom1_size, graph)
	
	#Second chromosome
	implement_second_chromosome_to_graph(unique_id, extracted, chrom2_size, graph)

	genome2_start_index_graph = length_common + length_unique_genome1

	#Set start index node in graph
	graph.set_start_index(genome1_start_index_graph, genome2_start_index_graph)

	print "\nGraph structure is finished\n"
from graph_v17 import *
from initialize import *

import time
import sys

def multiple_positions_from_files(graph, positions_filename_genome1, positions_filename_genome2, verbose):
	"""
	Initalize values from files and start algorithm

	Args:
		graph 						 (class): The graph class
		positions_filename_genome1  (string): Filename of positions for genome1
		positions_filename_genome2  (string): Filename of positions for genome2
		verbose 			   	   (boolean): Print out steps
	"""

	#Get a dictonary from the test data, where chromosome is the key
	dictionary_of_positions1 = parse_combined_file(positions_filename_genome1)
	dictionary_of_positions2 = parse_combined_file(positions_filename_genome2)
	print "Positions from file is finished read"

	#Max steps
	max_steps = eval(raw_input("Max steps: "))
	if (not isinstance(max_steps, int) or (max_steps < 0)):
		print "\nMax steps is not a positive integer"
		return

	count = 0

	#Time to see how long it took to run
	start_time = time.time()
	for chrom1, desired_positions in sorted(dictionary_of_positions1.items()):
		print "-----------Chromosome ", chrom1, "-----------"
		for position1 in desired_positions:
			print "Position genome1:", position1

			#Send the whole dictionary in the graph
			position2_dict = sorted(dictionary_of_positions2.items())
			result = graph.graph_based_liftover(position1, chrom1, position2_dict,
												max_steps, verbose)
			count += result
			print

			#return #Use max step 20 and check from 0 -> 2 -> 7 -> 1 (stop here)
		print

	print "*************"
	print ("max steps = %d" % max_steps) 
	print ("Time: %.2f minutes") % ((time.time() - start_time) / 60.0)
	print "Result:", count

def run(verbose):
	"""
	Run the real datasets

	Args:
		verbose (boolean): Print out steps
	"""

	#Read in 19 axt files
	filename = []
	for i in range(1,20):
		###CHANGE THIS TO FIND YOUR PATH
		filename.append(open("data/real_data/chr" + str(i) + ".mm10.hg19.net.axt"))
	
	###CHANGE THIS TO FIND YOUR PATH
	genome1_size = open("data/real_data/mm10_chrom_size.txt", "r")
	genome2_size = open("data/real_data/hg19_chrom_size.txt", "r")

	###CHANGE THIS TO FIND YOUR PATH
	positions_filename_genome1 = open("data/real_data/positions/mouse.bed")
	positions_filename_genome2 = open("data/real_data/positions/human.bed")
	
	#Initialize graph
	graph = Graph()

	#Make the graph
	start_time = time.time()
	initialize_file(filename, genome1_size, genome2_size, graph)
	print "*************Read the files in *****************"
	print ("Time: %.2f minutes") % ((time.time() - start_time) / 60.0)
	
	#print "Nodes in the graph:", len(graph.graph)

	multiple_positions_from_files(graph, positions_filename_genome1, 
								  positions_filename_genome2, verbose)

def test_run(verbose):
	"""
	Run the test datasets

	Args:
		verbose (boolean): Print out steps
	"""

	###CHANGE THIS TO FIND YOUR PATH
	filename1 = open("test_data/main_chr1.txt", "r")
	filename2 = open("test_data/main_chr2.txt", "r")
	filename = [filename1, filename2]
	
	###CHANGE THIS TO FIND YOUR PATH
	genome1_size = open("test_data/mouse_size.txt", "r")
	genome2_size = open("test_data/homo_size.txt", "r")
	
	###CHANGE THIS TO FIND YOUR PATH
	positions_filename_genome1 = open("test_data/mouse.txt")
	positions_filename_genome2 = open("test_data/human.txt")
	
	#Initialize graph
	graph = Graph()

	#Make the graph
	initialize_file(filename, genome1_size, genome2_size, graph)
	
	#print "Nodes in the graph:", len(graph.graph)

	multiple_positions_from_files(graph, positions_filename_genome1, positions_filename_genome2, verbose)

if __name__=="__main__":
	sys.setrecursionlimit(30000)
	
	#choose = raw_input("Verbose?: ")
	#verbose = False
	verbose = True

	#if (choose == "Y" or choose == "y"):
	#	verbose = True
	
	#run(verbose)
	test_run(verbose)
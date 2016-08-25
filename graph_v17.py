from constants import *

class Graph:
	def __init__(self):
		"""
		Initialize the graph
		"""

		self.graph = [] 			#Contain all the nodes. Parsing file do the work
		self.start_node = None 		#Define the start node id

		self.count = 0				#Add one if position where found
		
		self.position_dict = {}		#This variable contain the human genome, where 
									#the key is chr and value list of positions
		self.node_visited = []		#Nodes visited where node id is stored
		self.mode_forward = False	#Value to indicate forward method is commencing
		self.found = False			#If a position is found within max step, set to True

		self.node_front = {}		#Dictionary of shortest path from front of the node from the start node
		self.node_back = {}			#Dictionary of shortest path from the back of the node from the start node

	def add_common(self, node):
		"""
		Add common to the graph

		Args:
			node (class): Class node
		"""

		self.graph.append(node)

	def add_unique(self, list_node):
		"""
		Add the unique node to the graph

		Args:
			list_node (list): List of class nodes
		"""

		self.graph += list_node

	def set_start_index(self, genome1_start_index, genome2_start_index):
		"""
		Set the start index for genome1 and genome2 in the graph

		Args:
			genome1_start_index (int): Integer for indexation
			genome2_start_index (int): Integer for indexation
		"""

		self.genome1_start_index = genome1_start_index
		self.genome2_start_index = genome2_start_index

	def get_genome1_chromosome(self):
		"""
		Get the chromosome of genome1

		Returns:
			The chromosome of genome1
		"""

		return self.graph[self.genome1_start_index].chromosome_genome1

	def graph_based_liftover(self, position_genome1, chrom1, position_dict_read, max_steps, verbose):
		"""
		Initialize values before starting the algorithm

		Args:
			position_genome1       (int): Position to genome1
			chrom1 			    (string): Chromosome name for genome1
			position_dict_read 	  (list): Dictionary of chromosome with a list of position
			max_steps 		 	   (int): The max step compared to the current length
			verbose 		   (boolean): Print out steps

		Returns:
			The count value
		"""
		
		genome1_node = self.find_genome1_node(chrom1, position_genome1)

		if (genome1_node == None):
			return 0

		# if (verbose):
		# 	print "-------", genome1_node.id, "-------"

		self.position_dict = position_dict_read
		self.start_node = genome1_node.id
		
		#Check previous nodes in current node
		current_length = genome1_node.go_left_genome1(position_genome1)
		node_list = self.remove_start_or_end(genome1_node.previous)
		red_once = False
		direction = "previous"
		mode = "backward"
		previous_length = -1
		visit = []
		chrom2 = None

		#If recursion comes back to start node, no point going further. Let the forward do the job
		self.node_visited.append(genome1_node.id)

		self.find_match(current_length, genome1_node, node_list, chrom2, red_once, 
						max_steps, direction, mode, previous_length, visit, verbose)
		
		#Add the node in previous edge in current node
		self.node_backwards(genome1_node)
		
		#Check next nodes for current node
		current_length = genome1_node.go_right_genome1(position_genome1)
		node_list = self.remove_start_or_end(genome1_node.next)
		red_once = False
		direction = "next"
		mode = "forward"
		previous_length = 0
		self.mode_forward = True
		visit = []
		chrom2 = None

		#if (verbose):
		#	print "****************FORWARD****************"

		#self.node_visited.append(genome1_node.id)
	
		self.find_match(current_length, genome1_node, node_list, chrom2, red_once, 
						max_steps, direction, mode, previous_length, visit, verbose)
		
		return_value = self.count

		#Reset all variables
		self.reset_global_values()

		return return_value

	def find_match(self, current_length, current_node, node_list, chrom2, red_once, max_steps, 
				   direction, mode, previous_length, visit, verbose):
		"""
		This method start the algorithm

		Args:
			current_length 			(int): The current length so far
			current_node 		  (class): A node class
			node_list 			   (list): List of nodes (Can either be previous or next edge)
			chrom2 	 			 (string): Chromosome name if a red path has been visited
			red_once 			(boolean): True if a red path has been cross one time, else False
			max_steps 				(int): The max step compared to the current length
			direction 			 (string): The direction going through the nodes
			mode 				 (string): Backward or forward edge to use
			previous_length 		(int): The length before the on current node
			visit 				   (list): Visited nodes added during searching
			verbose 			(boolean): Print out steps

		Returns:
			None to end the recursive method
		"""
		
		# if (verbose):
		# 	print "-----------> Current ID:", current_node.id

		if (self.found):
			# if (verbose):
			# 	print "ALL FOUND, ENDING ALGORITHM"
			return
		
		#If current_length is lower than max_steps, keep going
		if (current_length <= max_steps):
			#Check if position is in node even the max step is not achived
			self.check_positions_list(current_node, verbose)

			if (self.found):
				return

			if (mode == "backward"):
				for prev_node in node_list:
					# if (verbose):
					# 	print "Update front on", prev_node.id, mode
					# 	print "MODE: Backward | prev:", current_node.id, "->" , prev_node.id , ":", current_length + 1 #UPDATE (front)

					if (not self.check_front(prev_node, current_length + 1)):
						# if (verbose):
						# 	print "MODE: Backward | prev (Check front is longer):", current_node.id
						pass
					else:
						temp_visit = visit + [current_node.id]
						self.find_match_backward(current_length, current_node, prev_node, chrom2, red_once, 
												 max_steps, direction, "previous", mode, previous_length, 
												 temp_visit, verbose)

					if (self.found):
						# if (verbose):
						# 	print "ALL self.found, ENDING:", current_node.id
						return

					for next_node in self.remove_start_or_end(current_node.next):
						if (current_node.id in self.node_visited):
							#This test is to assure the backward don't go forward from the start node. 
							#The calculation will be wrong in that case
							# if (verbose):
							# 	print "____________VISITED START NODE BACKWARD____________:", current_node.id, "| FORWARD TAKE CARE OF THIS:", next_node.id
							continue

						if (next_node.id in visit):
							# if (verbose):
							# 	print "____________VISITED BEFORE BACKWARD____________:", next_node.id
							continue

						#If a red path has been crossed, check if the next node has the same chromosome
						if (red_once and next_node.genome2 and next_node.chromosome_genome2 != chrom2):
							# if (verbose):
							# 	print "Current node:", current_node.id
							# 	print "Next node:", next_node.id, next_node.chromosome_genome2
							# 	print "YOU DON'T HAVE WHAT I NEED, BUT I HAVE TO GO THE METHOD ANYWAY"
							continue

						# if (verbose):
						# 	print "Update back on", next_node.id
						
						if (self.mode_forward):
							temp_length = current_length + 1
							# if (verbose):
							# 	print "Mode: Backward1 | next:", current_node.id, "->",  next_node.id, ":", temp_length #UPDATE (back)
						else:
							temp_length = current_length - previous_length + 1
							# if (verbose):
							# 	print "Mode: Backward2 | next:", current_node.id, "->",  next_node.id, ":", temp_length #UPDATE (back)

						if (not self.check_back(next_node, temp_length)):
							# if (verbose):
							# 	print "MODE: Backward | next (Check front is longer):", current_node.id
							continue

						temp_visit = visit + [current_node.id]
						self.find_match_forward(current_length, current_node, next_node, chrom2, 
												red_once, max_steps, direction, "next", mode, 
												previous_length, temp_visit, verbose)
			else:
				for next_node in node_list:
					# if (verbose):
					# 	print "Update back on", next_node.id, mode
					# 	print "Mode: Forward | next:", current_node.id, "->", next_node.id, ":", current_length + 1 #UPDATE (back)
					
					if (not self.check_back(next_node, current_length + 1)):
						# if (verbose):
						# 	print "MODE: Forward | next (Check front is longer):", current_node.id
						pass
					else:
						temp_visit = visit + [current_node.id]
						self.find_match_forward(current_length, current_node, next_node, chrom2,
												red_once, max_steps, direction, "next", mode,
												previous_length, temp_visit, verbose)
						
					if (self.found):
						# if (verbose):
						# 	print "ALL self.found, ENDING:", current_node.id
						return

					for prev_node in self.remove_start_or_end(current_node.previous):
						if (current_node.id in self.node_visited):
							#This test is to assure the backward don't go backward from the start node.
							#The calculation will be wrong in that case
							# if (verbose):
							# 	print "____________VISITED START NODE FORWARD____________:", current_node.id, "| BACKWARD ALREADY CHECKED THIS:", prev_node.id
							continue

						if (current_node.id in visit):
							# if (verbose):
							# 	print "____________VISITED B4444____________:", current_node.id
							continue

						if (red_once and prev_node.genome2 and prev_node.chromosome_genome2 != chrom2): #CHECK HERE!POPOFMODJOFJDOIFJ
							# if (verbose):
							# 	print "Current node:", current_node.id
							# 	print "Previous node:", prev_node.id, prev_node.chromosome_genome2
							# 	print "YOU DON'T HAVE WHAT I NEED, BUT I HAVE TO GO THE METHOD ANYWAY"
							continue

						if (prev_node.id in visit):
							# if (verbose):
							# 	print "VISITED BEFORE:", prev_node.id
							continue

						# if (verbose):
						# 	print "Update front on", prev_node.id
						
						if (self.mode_forward):
							temp_length = current_length + 1
							
#####						#Don't remember why I did this...
							if (self.mode_forward):
								temp_length -= previous_length

							# if (verbose):
							# 	print "Mode: Forward1 | prev:", current_node.id, "->",  prev_node.id, ":", temp_length #UPDATE (front)
						else:
							temp_length = current_length - previous_length + 1
							# if (verbose):
							# 	print "Mode: Forward2 | prev:", current_node.id, "->",  prev_node.id, ":", temp_length #UPDATE (front)

						if (not self.check_front(prev_node, temp_length)):
							# if (verbose):
							# 	print "MODE: Forward | prev (Check front is longer):", current_node.id
							continue

						temp_visit = visit + [current_node.id]
						self.find_match_backward(current_length, current_node, prev_node, chrom2, 
												 red_once, max_steps, direction, "previous", mode,
												 previous_length, temp_visit, verbose)
		else:
			#Current_length is larger than max_steps, so calculate the current node
			# if (verbose):
			# 	print "THE END, CALCULATING:", current_node.id
			# 	print "LENGTH (END):", current_length

			self.calculate(current_length, current_node, chrom2, red_once, max_steps, 
						   direction, mode, previous_length, verbose)

			if (self.found):
				# if (verbose):
				# 	print "ALL FOUND, ENDING:", current_node.id
				return

			if (mode == "backward"):
				# if (verbose):
				# 	print "AFTER CALCULATE:", mode

				for next_node in self.remove_start_or_end(current_node.next):
					if (next_node.id in visit):
						# if (verbose):
						# 	print "____________VISITED BEFORE check next____________:", next_node.id
						continue

					# if (verbose):
					# 	print "Try going in:", current_node.id, next_node.id

					self.check_node_after_max_step(current_length, current_node, next_node, chrom2,
						  						   red_once, max_steps, direction, mode, previous_length, 
					  							   verbose)
			else:
				# if (verbose):
				# 	print "AFTER CALCULATE:", mode

				for prev_node in self.remove_start_or_end(current_node.previous):
					if (prev_node.id in visit):
						# if (verbose):
						# 	print "____________VISITED BEFORE check prev____________:", prev_node.id
						continue

					# if (verbose):
					# 	print "Try going in:", current_node.id, prev_node.id,

					self.check_node_after_max_step(current_length, current_node, prev_node, chrom2,
						  						   red_once, max_steps, direction, mode, previous_length, 
					  							   verbose)

	def find_match_backward(self, current_length, current_node, prev_node, chrom2, red_once, max_steps, 
							current_direction, next_direction, mode, previous_length, visit, verbose):
		"""
		Use the previous edges on the current node

		Args:
			current_length 			(int): The current length so far
			current_node 		  (class): A node class
			prev_node 			  (class): The previous node of the current node
			chrom2 	 			 (string): Chromosome name if a red path has been visited
			red_once 			(boolean): True if a red path has been cross one time, else False
			max_steps 				(int): The max step compared to the current length
			current_direction 	 (string): The direction currently used now
			next_direction 		 (string): The direction on the next step
			mode 				 (string): Backward or forward edge to use
			previous_length 		(int): The length before the on current node
			visit 				   (list): Visited nodes added during searching
			verbose 			(boolean): Print out steps
		"""

		# if (verbose):
		# 	print "FIND MATCH BACKWARD"
		# 	print "RED visit before?:", red_once
		# 	print "Current Node:", current_node.id
		# 	print "Prev node:", prev_node.id
		
		if (current_node.path_type_backward(prev_node) == "blue_path"):
			# if (verbose):
			# 	print "BLUE (PREVIOUS)"
			# 	print "CURRENT LENGTH:", current_length

			if (mode == "forward"):
				#Going from forward to backward needs to reset the current length and add
				#the next length
				new_length = current_length - previous_length + prev_node.get_path_length_genome1() + 1
				mode = "backward"
				current_direction = "previous"
			else:
				previous_length = prev_node.get_path_length_genome1()
				new_length = current_length + previous_length + 1

				# if (verbose):
				# 	print "Added length:", previous_length

			# if (verbose):
			# 	print "NEW LENGTH (BACKWARD):", new_length
			# 	print "BACKWARD:", mode

			temp_list = self.remove_start_or_end(prev_node.previous)
			self.find_match(new_length, prev_node, temp_list, chrom2, red_once, max_steps, 
							current_direction, mode, previous_length, visit, verbose)
		else:
			# if (verbose):
			# 	print "RED (PREVIOUS)"
			# 	print "CURRENT LENGTH:", current_length

			if (red_once and current_node.chromosome_genome2 == chrom2):
				if (mode == "forward"):
					new_length = current_length - previous_length +  prev_node.get_path_length_genome2() + 1
					mode = "backward"
					current_direction = "previous"
				else:
					previous_length = prev_node.get_path_length_genome2()
					new_length = current_length + previous_length + 1
					# if (verbose):
					# 	print "Added length:", previous_length
				
				# if (verbose):
				# 	print "NEW LENGTH (BACKWARD):", new_length
				# 	print "BACKWARD:", mode

				temp_red_once = True
				temp_list = self.remove_start_or_end(prev_node.previous)
				self.find_match(new_length, prev_node, temp_list, chrom2, temp_red_once, 
								max_steps, current_direction, mode, previous_length, visit, verbose)
			else:
				if (mode == "forward"):
					new_length = current_length - previous_length + prev_node.get_path_length_genome2() + 1
					mode = "backward"
					current_direction = "previous"
				else:
					previous_length = prev_node.get_path_length_genome2()
					new_length = current_length + previous_length + 1
					# if (verbose):
					# 	print "Added length:", previous_length

				# if (verbose):
				# 	print "NEW LENGTH (BACKWARD):", new_length
				# 	print "BACKWARD:", mode

				temp_red_once = True
				temp_list = self.remove_start_or_end(prev_node.previous)
				self.find_match(new_length, prev_node, temp_list, prev_node.chromosome_genome2, 
								temp_red_once, max_steps, current_direction, mode, 
								previous_length, visit, verbose)

	def find_match_forward(self, current_length, current_node, next_node, chrom2, red_once, max_steps, 
						   current_direction, next_direction, mode, previous_length, visit, verbose):
		"""
		Use the next edges on the current node

		Args:
			current_length 			(int): The current length so far
			current_node 		  (class): A node class
			next_node 			  (class): The next node of the current node
			chrom2 	 			 (string): Chromosome name if a red path has been visited
			red_once 			(boolean): True if a red path has been cross one time, else False
			max_steps 				(int): The max step compared to the current length
			current_direction 	 (string): The direction currently used now
			next_direction 		 (string): The direction on the next step
			mode 				 (string): Backward or forward edge to use
			previous_length 		(int): The length before the on current node
			visit 				   (list): Visited nodes added during searching
			verbose 			(boolean): Print out steps
		"""

		# if (verbose):
		# 	print "FIND MATCH FORWARD"
		# 	print "RED visit before?:", red_once
		# 	print "Current Node:", current_node.id
		# 	print "Next node:", next_node.id

		if (current_node.path_type_forward(next_node) == "blue_path"):
			# if (verbose):
			# 	print "BLUE (NEXT)"
			# 	print "CURRENT LENGTH:", current_length

			if (mode == "backward"):
				#Going from backward to forward needs to reset the current length and add
				#the next length
				new_length = current_length - previous_length + next_node.get_path_length_genome1() + 1 
				mode = "forward"
				current_direction = "next"
			else:
				previous_length = next_node.get_path_length_genome1()
				new_length = current_length + previous_length + 1
				# if (verbose):
				# 	print "Added length:", previous_length

			# if (verbose):
			# 	print "NEW LENGTH (BACKWARD):", new_length
			# 	print "FORWARD:", mode

			temp_list = self.remove_start_or_end(next_node.next)
			self.find_match(new_length, next_node, temp_list, chrom2, red_once, max_steps, 
							current_direction, mode, previous_length, visit, verbose)
		else:
			if (red_once and current_node.chromosome_genome2 == chrom2):
				# if (verbose):
				# 	print "Current node:", current_node.id, next_node.id
				# 	print "RED BEFORE, BUT YOU STILL HAVE THE SAME CHROMOSOME. GOING FOR IT:", current_node.chromosome_genome2, next_node.chromosome_genome2
				
				if (mode == "backward"):
					#Going from backward to forward needs to reset the current length and add
					#the next length
					new_length = current_length - previous_length + next_node.get_path_length_genome2() + 1
					mode = "forward"
					current_direction = "next"
				else:
					previous_length = next_node.get_path_length_genome2()
					new_length = current_length + previous_length + 1
					# if (verbose):
					# 	print "Added length:", previous_length
				
				# if (verbose):
				# 	print "NEW LENGTH (BACKWARD):", new_length
				# 	print "FORWARD:", mode

				temp_red_once = True
				temp_list = self.remove_start_or_end(next_node.next)
				self.find_match(new_length, next_node, temp_list, chrom2, temp_red_once, 
								max_steps, current_direction, mode, previous_length, visit, verbose)
			else:
				# if (verbose):
				# 	print "RED (NEXT)"
				# 	print "CURRENT LENGTH:", current_length

				if (mode == "backward"):
					#Going from backward to forward needs to reset the current length and add
					#the next length
					new_length = current_length - previous_length + next_node.get_path_length_genome2() + 1
					mode = "forward"
					current_direction = "next"
				else:
					previous_length = next_node.get_path_length_genome2()
					new_length = current_length + previous_length + 1
					# if (verbose):
					# 	print "Added length:", previous_length
				
				# if (verbose):
				# 	print "NEW LENGTH (BACKWARD):", new_length
				# 	print "FORWARD:", mode

				temp_red_once = True
				temp_list = self.remove_start_or_end(next_node.next)
				self.find_match(new_length, next_node, temp_list, next_node.chromosome_genome2,
								temp_red_once, max_steps, current_direction, mode, 
								previous_length, visit, verbose)


	def check_positions_list(self, current_node, verbose):
		"""
		Check if the positions contain in the node

		Args:
			current_node  (class): The current node we want to check
			verbose 	(boolean): Print out steps
		"""

		if (not current_node.genome2):
			return

		positions_list = self.check_position(current_node)
		if (not positions_list):
			# if (verbose):
			# 	print "Can't find position in current node"
			return

		for pos in positions_list:
			if (current_node.check_position_genome2(pos)):
				self.count += 1
				self.found = True
				
				if (verbose):
					print "--> Found", current_node.chromosome_genome2, "in position", pos
				# if (verbose):
				# 	print "---------------------------------> FOUND (check): ", pos, "(Inside the node)"

	def calculate(self, current_length, current_node, chrom2, red_once, max_steps, 
				  direction, mode, previous_length, verbose):
		"""
		Calculate the remaning steps to check if positions is in the current node

		Args:
			current_length 			(int): The current length so far
			current_node 		  (class): The current node we want to check
			chrom2 	 			 (string): Chromosome name if a red path has been visited
			red_once 			(boolean): True if a red path has been cross one time, else False
			max_steps 				(int): The max step compared to the current length
			direction 	 		 (string): The direction currently used now
			mode 				 (string): Backward or forward edge to use
			previous_length 		(int): The length before the on current node
			verbose 			(boolean): Print out steps
		"""

		# if (verbose):
		# 	print "IN CALC:", current_node.id, current_length

		#If current node don't contain genome2, no need to further investigate
		if (not current_node.genome2):
			return

		#Subtract the current length by the current nodes length, because of the recursive
		current_length -= current_node.get_path_length_genome2()
		
		positions_list = self.check_position(current_node)
		if (not positions_list):
			# if (verbose):
			# 	print "Can't find position in current node"
			return

		for pos in positions_list:
			if (not current_node.check_position_genome2(pos)):
				continue	

			if (direction == "previous"):
				rest_length = current_node.go_left_genome2(pos)
			else:
				rest_length = current_node.go_right_genome2(pos)

			new_length = current_length + rest_length
			
			if (new_length <= max_steps):
				self.count += 1
				self.found = True
				
				if (verbose):
					print "Found", current_node.chromosome_genome2, "in position", pos
				# if (verbose):
				# 	print "---------------------------------> FOUND (check): ", pos, new_length


	def check_node_after_max_step(self, current_length, current_node, prev_next_node, chrom2, red_once, 
								  max_steps, direction, mode, previous_length, verbose):
		"""
		Check on the current node the opposite direction if there is possible path. Assume there will
		maximum be one node to check

		Args:
			current_length 			(int): The current length so far
			current_node 		  (class): A node class
			prev_next_node 		  (class): Can either be a previous or next node
			chrom2 	 			 (string): Chromosome name if a red path has been visited
			red_once 			(boolean): True if a red path has been cross one time, else False
			max_steps 				(int): The max step compared to the current length
			direction 	 		 (string): The direction currently used now
			mode 				 (string): Backward or forward edge to use
			previous_length 		(int): The length before the on current node
			verbose 			(boolean): Print out steps
		"""

		# if (verbose):
		# 	print "METHOD CHECK, ", mode, direction

	 	current_length -= previous_length
	 	
 		# if (verbose):
 		# 	print "CHECK CALCULATE", current_node.id, prev_next_node.id

		if (not prev_next_node.genome2):
			# if (verbose):
			# 	print "In check: Don't contain genome2:", prev_next_node.id
			return

		if (red_once and prev_next_node.chromosome_genome2 != chrom2):
			# if (verbose):
			# 	print "Current node:", current_node.id
			# 	print "Checking next node:", prev_next_node.id, prev_next_node.chromosome_genome2
			# 	print "I'M IN CHECK METHOD, LAST STEP. FOUND CHR DIFFERENT. ENDING CALCULATE"
			return

		# if (verbose):
		# 	print "CHECK CALC (NO MORE STEP)-->:", current_node.id, "->", prev_next_node.id, ":", current_length + 1 #UPDATE (depends)

		if (mode == "forward"):
			# if (verbose):
			# 	print "Update front on ", prev_next_node.id

			direction = "previous"

			if (not self.check_front(prev_next_node, current_length + 1)):
				# if (verbose):
				# 	print "MODE: CALCULATE | FRONT (Check front is longer):", prev_next_node.id
				return
		else:
			# if (verbose):
			# 	print "Update back on ", prev_next_node.id

			direction = "next"

			if (not self.check_back(prev_next_node, current_length + 1)):
				# if (verbose):
				# 	print "MODE: CALCULATE | BACK (Check back is longer):", prev_next_node.id
				return

		temp_length = prev_next_node.get_path_length_genome2() + 1
		current_length += temp_length

		# if (verbose):
		# 	print "FOUND NODE:", prev_next_node.id

		self.calculate(current_length, prev_next_node, chrom2, red_once, max_steps, direction, 
					   mode, prev_next_node.get_path_length_genome2(), verbose)
	
	def check_front(self, current_node, step):
		"""
		Check if the current node has been added before. If not, add the node with the current path
		length append as the value. Else check if the new path is shorter than the previous path.
		This method add to the front of the current node

		Args:
			current_node (class): A node class
			step 		   (int): The step taken to the current node

		Returns:
			True if path has not been visited before or new path is shorter than the previous one
			Else false
		"""

		if (current_node.id not in self.node_front):
			self.node_front[current_node.id] = step
			return True
		elif (step < self.node_front[current_node.id]):
			self.node_front[current_node.id] = step
			return True

		return False

	def check_back(self, current_node, step):
		"""
		Check if the current node has been added before. If not, add the node with the current path
		length append as the value. Else check if the new path is shorter than the previous path.
		This method add to the back of the current node

		Args:
			current_node (class): A node class
			step 		   (int): The step taken to the current node

		Returns:
			True if path has not been visited before or new path is shorter than the previous one
			Else false
		"""

		if (current_node.id not in self.node_back):
			self.node_back[current_node.id] = step
			return True
		elif (step < self.node_back[current_node.id]):
			self.node_back[current_node.id] = step
			return True

		return False

	def check_position(self, current_node):
		"""
		Check if the current_node contain the chromosome number

		Args:
			current_node (class): A node class

		Returns:
			True returns a list of the position
			False returns False
		"""

		for chrom, positions_list in self.position_dict:
			if (current_node.chromosome_genome2 == chrom):
				return positions_list

		return False

	def find_genome1_node(self, chrom1_genome1, position_genome1):
		"""
		Find the node to genome1 by the position_genome1

		Args:
			chrom1_genome1 	 (string): Chromosome name
			position_genome1 	(int): Position to genome1

		Returns:
			Return the class node. If not found, then return None
		"""

		for node in self.graph:
			if (node.genome1 and node.chromosome_genome1 == chrom1_genome1):
				if (node.check_position_genome1(position_genome1)):
					return node

		return None

	def remove_start_or_end(self, node_list):
		"""
		Remove the start and end node in the edges list

		Args:
			node_list (list): List of nodes

		Returns:
			A new list without start and end node
		"""

		if ("start" in node_list):
			node_list.remove("start")

		if ("end" in node_list):
			node_list.remove("end")

		return node_list
	
	def node_backwards(self, current_node):
		for node in self.remove_start_or_end(current_node.previous):
			if (node.id not in self.node_visited):
				self.node_visited.append(node.id)

	def reset_global_values(self):
		"""
		Reset count, positions_list and node_visited
		"""

		self.count = 0
		self.node_visited = []
		self.mode_forward = False
		self.node_front = {}
		self.node_back = {}
		self.found = False
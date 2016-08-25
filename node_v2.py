class Node:
	def  __init__(self, id):
		"""
		Initialize values

		Args:
			id (int): Identity to the node
		"""

		self.id = id

		self.previous = []
		self.next = []

		self.path_stop = []

		self.genome1 = False
		self.genome2 = False

	#-------------------------------SET GENOME METHOD-------------------------------
	def set_genome1(self, chromosome, start, end):
		"""
		Set the parameters for genome1 (main chromosome)

		Args:
			chromosome (string): Chromosome name with integer (ex. chr1)
			start 		  (int): Start position
			end 		  (int): End position
		"""

		self.chromosome_genome1 = chromosome
		self.start_genome1 = start
		self.end_genome1 = end
		self.genome1 = True

	def set_genome2(self, chromosome, start, end):
		"""
		Set the parameters for genome2 (aligning chromosome)

		Args:
			chromosome (string): Chromosome name with integer (ex. chr1)
			start 		  (int): Start position
			end 		  (int): End position
		"""
		
		self.chromosome_genome2 = chromosome
		self.start_genome2 = start
		self.end_genome2 = end
		self.genome2 = True

	#-------------------------------SET EDGE POINTER-------------------------------
	def set_previous_pointer(self, previous_pointer):
		"""
		Set previous node

		Args:
			previous_pointer (class): Node class
		"""
		
		self.previous.append(previous_pointer)

	def set_next_pointer(self, next_pointer):
		"""
		Set next node

		Args:
			next_pointer (class): Node class
		"""
		
		self.next.append(next_pointer)

	#-------------------------------GET NODE METHOD-------------------------------
	def get_next_genome1(self):
		"""
		Get the next node of genome1 (main chromosome)

		Returns:
			The next node class
		"""

		return self.next[0]

	def get_next_genome2(self):
		"""
		#THIS METHOD IS NOT USED AT ALL
		Get the next node of genome2 (aligning chromosome)

		Returns:
			The next node class
		"""

		index = 1
		if (len(self.next) != 2):
			index = 0

	 	return self.next[index]

	#-------------------------------GET LENGTH METHOD-------------------------------
	def get_path_length_genome1(self):
		"""
		Get the length of genome1

		Returns:
			Length of genome1
		"""
		
		return self.end_genome1 - self.start_genome1

	def get_path_length_genome2(self):
		"""
		Get the length of genome2

		Returns:
			Length of genome2
		"""

		return self.end_genome2 - self.start_genome2

	#----------------------------CALCULATE STEPS GENOME1----------------------------
	def go_left_genome1(self, position):
		"""
		Calculate the steps to the left for genome1 (main chromosome)

		Args:
			position (int): Integer value

		Returns:
			Steps from position to start position
		"""

		return position - self.start_genome1

	def go_right_genome1(self, position):
		"""
		Calculate the steps to the right for genome1 (main chromosome)

		Args:
			position (int): Integer value

		Returns:
			Steps from position to end position
		"""

		return self.end_genome1 - position

	#----------------------------CALCULATE STEPS GENOME2----------------------------
	def go_right_genome2(self, position):
		"""
		Calculate the steps from start position of the node to the position

		Args:
			position (int): Integer value

		Returns:
			Steps from position to start position
		"""

		return position - self.start_genome2

	def go_left_genome2(self, position):
		"""
		Calculate the steps from end position of the node to the position

		Args:
			position (int): Integer value

		Returns:
			Steps from position to end position
		"""

		return self.end_genome2 - position

	#--------------------------------CHECK POSITION--------------------------------
	def check_position_genome1(self, position):
		"""
		Check if position exists in this node for genome1 (main chromosome)

		Args:
			position (int): Integer value

		Returns:
			True -> Position is in this node
			False -> Else
		"""

		if ((self.start_genome1 <= position) and (position <= self.end_genome1)):
			return True
		return False

	def check_position_genome2(self, position):
		"""
		Check if position exists in this node for genome2 (aligning chromosome)

		Args:
			position (int): Integer value

		Returns:
			True -> Position is in this node
			False -> Else
		"""

		if ((self.start_genome2 <= position) and (position <= self.end_genome2)):
			return True
		return False
	
	#---------------------------------PATH ILLEGAL---------------------------------
	def path_not_allowed(self, node_id):
		"""
		#THIS METHOD IS NOT USED AT ALL
		Add the node id of illegal path

		Args:
			node_id (int): Node id
		"""

		self.path_stop.append(node_id)

	#----------------------------------PATH TYPE----------------------------------
	def path_type_backward(self, node):
		"""
		Identify the path that can be taken. This method is used by find_match_backward

		Args:
			Node (class): The node class

		Returns:
			Blue or red path
		"""
		
		if (self.genome1 and self.genome2 and node.genome1 and node.genome2):
			if (self.start_genome1 - node.end_genome1 == 1):
				return "blue_path"
			else:
				return "red_path"
		elif (self.genome1 and node.genome1):
			return "blue_path"
		elif (self.genome2 and node.genome2):
			return "red_path"
		elif (node.genome1):
			return "blue_path"
		elif (node.genome2):
			return "red_path"

	def path_type_forward(self, node):
		"""
		Identify the path that can be taken. This method is used by find_match_forward

		Args:
			Node (class): The node class

		Returns:
			Blue or red path
		"""
		
		if (self.genome1 and self.genome2 and node.genome1 and node.genome2):
			if (node.start_genome1 - self.end_genome1 == 1):
				return "blue_path"
			else:
				return "red_path"
		elif (self.genome1 and node.genome1):
			return "blue_path"
		elif (self.genome2 and node.genome2):
			return "red_path"
		elif (node.genome1):
			return "blue_path"
		elif (node.genome2):
			return "red_path"

	#-----------------------------CALCULATE INNER STEP-----------------------------
	def calculate_inner_step(self, position_genome1, position_genome2):
		"""
		#THIS METHOD IS NEVER USED
		Calculate the step from position in genome1 and genome2, when they have different 
		length. This method is only called when position is in the same node

		Args:
			position_genome1 (int): Position to genome1
			position_genome2 (int): Position to genome2

		Returns:
			The difference in length
		"""

		length_genome1 = abs(position_genome1 - self.start_genome1)
		length_genome2 = abs(position_genome2 - self.start_genome2)

		return abs(length_genome1 - length_genome2)
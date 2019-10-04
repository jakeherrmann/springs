
from .node      import Node
from .spring    import Spring
from .boundary  import Boundary
from .springFileIO import FileVariable, FileFormat
from pathlib import Path

import operator
import struct
import os

class SpringNetwork:
	relops = {
		'>' : operator.gt ,
		'>=': operator.ge ,
		'<' : operator.lt ,
		'<=': operator.le ,
		'==': operator.eq ,
		'!=': operator.ne }

	def __init__(self, num_dimensions=2, precision='float'):
		self.num_nodes = 0
		self.num_springs = 0
		self.num_boundaries = 0
		self.num_dimensions = num_dimensions
		self.precision = precision
		self.num_stiffness_tension = 1
		self.num_stiffness_compression = 0
		self.nodes = []
		self.springs = []
		self.boundaries = []
		self.dir_solver_input = None
		self.dir_solver_output = None

	def setup(self, nodes=None, springs=None, boundaries=None):
		if nodes is not None:
			self.nodes = nodes
			self.num_nodes = len(nodes)
		if springs is not None:
			self.springs = springs
			self.num_springs = len(springs)
			self.connect_springs_nodes()
		if boundaries is not None:
			self.boundaries = boundaries
			self.num_boundaries = len(boundaries)
			self.connect_boundaries_nodes()

	def solve(self):
		self.write_spring_network(self.dir_solver_input, use_solver_format=True)
		self.run_spring_solver(self.dir_solver_input, self.dir_solver_output)
		self.read_spring_network(self.dir_solver_output, use_solver_format=True)

	def run_spring_solver(self, dir_input=None, dir_output=None):
		if dir_input is None:
			dir_input = self.dir_solver_input
		exe_springs_solver = Path('.') / 'springs_solver' / 'springs_solver.exe'
		sys_command = str(exe_springs_solver) + ' ' + str(dir_input)
		if dir_output is not None:
			sys_command += ' ' + str(dir_output)
		os.system(sys_command)

	def write_spring_network(self, dir_input=None, use_solver_format=False):
		if dir_input is None:
			dir_input = self.dir_solver_input
		if not dir_input.exists():
			dir_input.mkdir(parents=True)
		write_nodes = self.nodes
		write_springs = self.springs if not use_solver_format else [ spring for spring in self.springs if not spring.broken ]
		file_name = dir_input / 'network_parameters.txt'
		with open(file_name,'wt') as file:
			file.write('{:d}\n'.format(len(write_nodes)))
			file.write('{:d}\n'.format(len(write_springs)))
			file.write('{:s}\n'.format(self.precision))
			file.write('{:d}\n'.format(self.num_dimensions))
			file.write('{:d}\n'.format(self.num_stiffness_tension))
			file.write('{:d}\n'.format(self.num_stiffness_compression))
		file_name = dir_input / 'network_nodes.dat'
		file_format = Node.get_file_format(
			self.num_dimensions,
			self.precision,
			use_solver_format=use_solver_format)
		file_format.write_binary_file(file_name, write_nodes)
		file_name = dir_input / 'network_springs.dat'
		file_format = Spring.get_file_format(
			self.precision,
			self.num_stiffness_tension,
			self.num_stiffness_compression,
			use_solver_format=use_solver_format)
		file_format.write_binary_file(file_name, write_springs)
		if not use_solver_format:
			file_name = dir_input / 'network_boundaries.dat'
			file_format = Boundary.get_file_format(
				self.num_dimensions)
			file_format.write_binary_file(file_name, self.boundaries)

	def read_spring_network(self, dir_output=None, use_solver_format=False, reinitialize=False):
		if dir_output is None:
			dir_output = self.dir_solver_output
		file_name = dir_output / 'network_parameters.txt'
		if reinitialize:
			with open(file_name,'rt') as file:
				self.num_nodes                 = int(file.readline())
				self.num_springs               = int(file.readline())
				self.precision                 = file.readline().rstrip('\n')
				self.num_dimensions            = int(file.readline())
				self.num_stiffness_tension     = int(file.readline())
				self.num_stiffness_compression = int(file.readline())
			self.nodes   = [ Node(self.num_dimensions) for count in range(self.num_nodes) ]
			self.springs = [ Spring() for count in range(self.num_springs) ]
			if not use_solver_format:
				self.boundaries = [ Boundary() for count in range(self.num_boundaries) ]
		read_nodes = self.nodes
		read_springs = self.springs if not use_solver_format else [ spring for spring in self.springs if not spring.broken ]
		file_name = dir_output / 'network_nodes.dat'
		file_format = Node.get_file_format(
			self.num_dimensions,
			self.precision,
			use_solver_format=use_solver_format)
		file_format.read_binary_file(file_name, read_nodes)
		file_name = dir_output / 'network_springs.dat'
		file_format = Spring.get_file_format(
			self.precision,
			self.num_stiffness_tension,
			self.num_stiffness_compression,
			use_solver_format=use_solver_format)
		file_format.read_binary_file(file_name, read_springs)
		if not use_solver_format:
			file_name = dir_output / 'network_boundaries.dat'
			file_format = Boundary.get_file_format(
				self.num_dimensions)
			file_format.read_binary_file(file_name, self.boundaries)
			for b in self.boundaries:
				b.force_directions = [ b.force_directions[i:i+self.num_dimensions] for i in range(0, len(b.force_directions), self.num_dimensions) ]
				b.displacements    = [ b.displacements   [i:i+self.num_dimensions] for i in range(0, len(b.displacements   ), self.num_dimensions) ]
		if reinitialize:
			self.connect_springs_nodes()
			self.connect_boundaries_nodes()

	def connect_springs_nodes(self):
		for spring in self.springs:
			spring.node_start_pointer = self.nodes[spring.node_start]
			spring.node_end_pointer   = self.nodes[spring.node_end  ]

	def connect_boundaries_nodes(self):
		if self.boundaries is not None:
			for boundary in self.boundaries:
				boundary.nodes_pointer = [ self.nodes[n] for n in boundary.nodes ]

	def reset_boundary_conditions(self, which_conditions='all'):
		if self.boundaries is not None:
			for boundary in self.boundaries:
				boundary.reset_conditions(which_conditions)

	def apply_boundary_conditions(self):
		for node in self.nodes:
			node.force = [0.0] * self.num_dimensions
		if self.boundaries is not None:
			for boundary in self.boundaries:
				boundary.apply_conditions()

	def apply_stretch(self, stretch, dimensions='all', boundary_indexes='all'):
		self.reset_boundary_conditions('displacements')
		if dimensions is not None:
			if dimensions=='all':
				dimensions = range(self.num_dimensions)
			if boundary_indexes is None:
				for dim in dimensions:
					self.apply_stretch_everywhere(stretch, dim)
			else:
				if boundary_indexes=='all':
					boundary_indexes = range(len(self.boundaries))
				for bi in boundary_indexes:
					boundary = self.boundaries[bi]
					for n in boundary.nodes_pointer:
						boundary.displacements.append( [ (stretch-1.0)*x if d in dimensions else 0.0 for d, x in enumerate(n.position) ] )
		self.apply_boundary_conditions()

	def apply_stretch_everywhere(self, stretch, dim):
		for node in self.nodes:
			node.position[dim] *= stretch

	def calc_spring_force(self):
		for spring in self.springs:
			spring.calc_force()

	def break_spring(self, break_variable, threshold, relop='>='):
		if hasattr(self.springs[0], break_variable):
			# now, set all nodes to not referenced
			# later, set the nodes of non-broken springs to referenced
			for node in self.nodes:
				node.referenced = False
			for spring in self.springs:
				if not spring.broken:
					value = getattr(spring, break_variable)
					value = value[0] if isinstance(value, list) else value
					spring.broken = SpringNetwork.relops[relop](value, threshold)
					if not spring.broken:
						spring.node_start_pointer.referenced = True
						spring.node_end_pointer.referenced = True

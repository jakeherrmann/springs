
from .node     import Node
from .spring   import Spring
from .boundary import Boundary
from .springFileIO import FileVariable, FileFormat
from pathlib import Path
import struct
import os

class SpringNetwork:
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
		self.dir_input = None
		self.dir_output = None

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
		self.write_spring_network(self.dir_input)
		self.run_spring_solver(self.dir_input, self.dir_output)
		self.read_spring_network(self.dir_output)

	def run_spring_solver(self, dir_input=None, dir_output=None):
		if dir_input is None:
			dir_input = self.dir_input
		exe_springs_solver = Path('.') / 'springs_solver' / 'springs_solver.exe'
		sys_command = str(exe_springs_solver) + ' ' + str(dir_input)
		if dir_output is not None:
			sys_command += ' ' + str(dir_output)
		os.system(sys_command)

	def write_spring_network(self, dir_input=None):
		if dir_input is None:
			dir_input = self.dir_input
		if not dir_input.exists():
			dir_input.mkdir(parents=True)
		file_name = dir_input / 'network_parameters.txt'
		with open(file_name,'wt') as file:
			file.write('{:d}\n'.format(self.num_nodes))
			file.write('{:d}\n'.format(self.num_springs))
			file.write('{:s}\n'.format(self.precision))
			file.write('{:d}\n'.format(self.num_dimensions))
			file.write('{:d}\n'.format(self.num_stiffness_tension))
			file.write('{:d}\n'.format(self.num_stiffness_compression))
		file_name = dir_input / 'network_nodes.dat'
		file_format = Node.get_file_format(self.num_dimensions, self.precision)
		file_format.write_binary_file(file_name, self.nodes)
		file_name = dir_input / 'network_springs.dat'
		file_format = Spring.get_file_format(self.precision, self.num_stiffness_tension, self.num_stiffness_compression)
		file_format.write_binary_file(file_name, self.springs)

	def read_spring_network(self, dir_output=None):
		if dir_output is None:
			dir_output = self.dir_output
		file_name = dir_output / 'network_parameters.txt'
		with open(file_name,'rt') as file:
			self.num_nodes                 = int(file.readline())
			self.num_springs               = int(file.readline())
			self.precision                 = file.readline().rstrip('\n')
			self.num_dimensions            = int(file.readline())
			self.num_stiffness_tension     = int(file.readline())
			self.num_stiffness_compression = int(file.readline())
		self.nodes   = [ Node(self.num_dimensions) for count in range(self.num_nodes  ) ]
		self.springs = [ Spring() for count in range(self.num_springs) ]
		file_name = dir_output / 'network_nodes.dat'
		file_format = Node.get_file_format(self.num_dimensions, self.precision)
		file_format.read_binary_file(file_name, self.nodes)
		file_name = dir_output / 'network_springs.dat'
		file_format = Spring.get_file_format(self.precision, self.num_stiffness_tension, self.num_stiffness_compression)
		file_format.read_binary_file(file_name, self.springs)
		self.connect_springs_nodes()

	def connect_springs_nodes(self):
		for spring in self.springs:
			spring.node_start_pointer = self.nodes[spring.node_start]
			spring.node_end_pointer   = self.nodes[spring.node_end  ]

	def connect_boundaries_nodes(self):
		for boundary in self.boundaries:
			boundary.nodes_pointer = [ self.nodes[n] for n in boundary.nodes ]

	def apply_boundary_conditions(self):
		for boundary in self.boundaries:
			boundary.apply_conditions()

	def apply_stretch(self, stretch, dim):
		for node in self.nodes:
			node.position[dim] *= stretch

	def calc_spring_force(self):
		for spring in self.springs:
			spring.calc_force()

	def break_spring_force(self, threshold):
		self.springs = [ spring for spring in self.springs if spring.force < threshold ]
		self.num_springs = len(self.springs)

	def break_spring_strain(self, threshold):
		self.springs = [ spring for spring in self.springs if spring.strain < threshold ]
		self.num_springs = len(self.springs)



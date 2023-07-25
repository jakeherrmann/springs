
from .node           import Node
from .spring         import Spring
from .structure      import Structure
from .structureGroup import StructureGroup
from .boundary       import Boundary
from .springFileIO   import FileVariable, FileFormat
from pathlib import Path

import springs_setup

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
		self.num_structures = 0
		self.num_structure_groups = 0
		self.num_boundaries = 0
		self.num_dimensions = num_dimensions
		self.precision = precision
		self.nodes = []
		self.springs = []
		self.structures = []
		self.structure_groups = []
		self.boundaries = []
		#
		self.dir_solver_input = None
		self.dir_solver_output = None
		self.solver_algorithm = 'newton'
		self.solver_objective = 'energy'
		self.solver_num_threads = 4
		self.solver_num_iter_save = 0
		self.solver_num_iter_print = 0
		self.solver_num_iter_max = 0
		self.solver_use_numerical_hessian = False
		self.solver_tolerance_change_objective = 1.0e-12
		self.solver_tolerance_sum_net_force = 1.0e-12
		self.solver_verbose = 1
		self.solver_parallel = 1
		#

	def setup(self, nodes=None, springs=None, structures=None, structure_groups=None, boundaries=None):
		if nodes is not None:
			self.nodes = nodes
			self.num_nodes = len(nodes)
		if springs is not None:
			self.springs = springs
			self.num_springs = len(springs)
			self.connect_springs_nodes()
		if structures is not None:
			self.structures = structures
			self.num_structures = len(structures)
			self.connect_structures_nodes_springs()
		if structure_groups is not None:
			self.structure_groups = structure_groups
			self.num_structure_groups = len(structure_groups)
			self.connect_structure_groups_structures()
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
		if dir_output is None:
			dir_output = self.dir_solver_output
		if not dir_output.exists():
			dir_output.mkdir(parents=True)
		exe_springs_solver = Path('.') / 'springs_solver' / 'springs_solver.exe'
		if not exe_springs_solver.exists():
			print('Could not find compiled solver executable.')
			print('Attempting to compile now.')
			springs_setup.main('')
		sys_command = str(exe_springs_solver) + ' ' + str(dir_input) + ' ' + str(dir_output)
		sys_command += ' --verbose ' + str(self.solver_verbose)
		sys_command += ' --parallel ' + str(self.solver_parallel)
		if self.solver_verbose:
			print(sys_command)
		os.environ['OMP_STACKSIZE'  ] = '64M'    ;
		os.environ['OMP_PLACES'     ] = 'cores'  ;
		os.environ['OMP_WAIT_POLICY'] = 'active' ;
		os.environ['OMP_PROC_BIND'  ] = 'true'   ;
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
			file.write(                 'num_points {:d}\n'.format(len(write_nodes)))
			file.write(                'num_springs {:d}\n'.format(len(write_springs)))
			file.write(                  'precision {:s}\n'.format(self.precision))
			file.write(             'num_dimensions {:d}\n'.format(self.num_dimensions))
			file.write(                  'algorithm {:s}\n'.format(self.solver_algorithm))
			file.write(                  'objective {:s}\n'.format(self.solver_objective))
			file.write(                'num_threads {:d}\n'.format(self.solver_num_threads))
			file.write(              'num_iter_save {:d}\n'.format(self.solver_num_iter_save))
			file.write(             'num_iter_print {:d}\n'.format(self.solver_num_iter_print))
			file.write(               'num_iter_max {:d}\n'.format(self.solver_num_iter_max))
			file.write(      'use_numerical_hessian {:d}\n'.format(self.solver_use_numerical_hessian))
			file.write( 'tolerance_change_objective {:e}\n'.format(self.solver_tolerance_change_objective))
			file.write(    'tolerance_sum_net_force {:e}\n'.format(self.solver_tolerance_sum_net_force))
		file_name = dir_input / 'network_nodes.dat'
		file_format = Node.get_file_format(
			self.num_dimensions,
			self.precision,
			use_solver_format=use_solver_format)
		file_format.write_binary_file(file_name, write_nodes)
		file_name = dir_input / 'network_springs.dat'
		file_format = Spring.get_file_format(
			self.precision,
			use_solver_format=use_solver_format)
		file_format.write_binary_file(file_name, write_springs)
		if not use_solver_format:
			file_name = dir_input / 'network_structures.dat'
			file_format = Structure.get_file_format()
			file_format.write_binary_file(file_name, self.structures)
			file_name = dir_input / 'network_structure_groups.dat'
			file_format = StructureGroup.get_file_format()
			file_format.write_binary_file(file_name, self.structure_groups)
			file_name = dir_input / 'network_boundaries.dat'
			file_format = Boundary.get_file_format(
				self.num_dimensions)
			file_format.write_binary_file(file_name, self.boundaries)

	def read_spring_network(self, dir_output=None, use_solver_format=False, reinitialize=False):
		if dir_output is None:
			dir_output = self.dir_solver_output
		if reinitialize:
			file_name = dir_output / 'network_parameters.txt'
			with open(file_name,'rt') as file:
				for line in file.readlines():
					arg = line.rstrip('\n').split()
					if   arg[0]=="num_points"                 : self.num_nodes                         =   int(arg[1])
					elif arg[0]=="num_springs"                : self.num_springs                       =   int(arg[1])
					elif arg[0]=="precision"                  : self.precision                         =       arg[1]
					elif arg[0]=="num_dimensions"             : self.num_dimensions                    =   int(arg[1])
					elif arg[0]=="algorithm"                  : self.solver_algorithm                  =       arg[1]
					elif arg[0]=="objective"                  : self.solver_objective                  =       arg[1]
					elif arg[0]=="num_threads"                : self.solver_num_threads                =   int(arg[1])
					elif arg[0]=="num_iter_save"              : self.solver_num_iter_save              =   int(arg[1])
					elif arg[0]=="num_iter_print"             : self.solver_num_iter_print             =   int(arg[1])
					elif arg[0]=="num_iter_max"               : self.solver_num_iter_max               =   int(arg[1])
					elif arg[0]=="use_numerical_hessian"      : self.solver_use_numerical_hessian      =  bool(arg[1])
					elif arg[0]=="tolerance_change_objective" : self.solver_tolerance_change_objective = float(arg[1])
					elif arg[0]=="tolerance_sum_net_force"    : self.solver_tolerance_sum_net_force    = float(arg[1])
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
			use_solver_format=use_solver_format)
		file_format.read_binary_file(file_name, read_springs)
		if not use_solver_format:
			file_name = dir_output / 'network_structures.dat'
			file_format = Structure.get_file_format()
			file_format.read_binary_file(file_name, self.structures)
			file_name = dir_output / 'network_structure_groups.dat'
			file_format = StructureGroup.get_file_format()
			file_format.read_binary_file(file_name, self.structure_groups)
			file_name = dir_output / 'network_boundaries.dat'
			file_format = Boundary.get_file_format(
				self.num_dimensions)
			file_format.read_binary_file(file_name, self.boundaries)
			for b in self.boundaries:
				b.force_directions = [ b.force_directions[i:i+self.num_dimensions] for i in range(0, len(b.force_directions), self.num_dimensions) ]
				b.displacements    = [ b.displacements   [i:i+self.num_dimensions] for i in range(0, len(b.displacements   ), self.num_dimensions) ]
		if reinitialize:
			self.connect_springs_nodes()
			self.connect_structures_nodes_springs()
			self.connect_structure_groups_structures()
			self.connect_boundaries_nodes()

	def connect_springs_nodes(self):
		for spring in self.springs:
			spring.node_start = self.nodes[spring.node_start_index]
			spring.node_end   = self.nodes[spring.node_end_index  ]
		for spring in self.springs:
			spring.adjacent_nodes.clear()
			spring.adjacent_springs.clear()
		for node in self.nodes:
			node.adjacent_nodes.clear()
			node.adjacent_springs.clear()
		for spring in self.springs:
			if not spring.broken:
				spring.node_start.adjacent_springs.append( spring )
				spring.node_start.adjacent_nodes.append( spring.node_end )
				spring.node_end.adjacent_springs.append( spring )
				spring.node_end.adjacent_nodes.append( spring.node_start )
		for spring in self.springs:
			if not spring.broken:
				for adj_spring in spring.node_start.adjacent_springs:
					if adj_spring is not spring and not adj_spring.broken:
						adj_spring.adjacent_springs.append( spring )
						adj_spring.adjacent_nodes.append( spring.node_end )
				for adj_spring in spring.node_end.adjacent_springs:
					if adj_spring is not spring and not adj_spring.broken:
						adj_spring.adjacent_springs.append( spring )
						adj_spring.adjacent_nodes.append( spring.node_start )

	def connect_structures_nodes_springs(self):
		for structure in self.structures:
			structure.nodes   = [ self.nodes[n]   for n in structure.nodes_indexes   ]
			structure.springs = [ self.springs[s] for s in structure.springs_indexes ]
		for node in self.nodes:
			node.structures.clear()
		for spring in self.springs:
			spring.structures.clear()
		for structure in self.structures:
			for node in structure.nodes:
				node.structures.append( structure )
			for spring in structure.springs:
				spring.structures.append( structure )
		for structure in self.structures:
			structure.adjacent_structures.clear()
		for structure in self.structures:
			for spring in structure.springs:
				for adj_structure in spring.structures:
					if adj_structure is not structure:
						structure.adjacent_structures.append( adj_structure )

	def connect_structure_groups_structures(self):
		if self.structure_groups is not None:
			for structure_group in self.structure_groups:
				structure_group.structures = [ self.structures[s] for s in structure_group.structures_indexes ]

	def connect_boundaries_nodes(self):
		if self.boundaries is not None:
			for boundary in self.boundaries:
				boundary.nodes = [ self.nodes[n] for n in boundary.nodes_indexes ]

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
					for n in boundary.nodes:
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
			any_new_breaks = False
			for node in self.nodes:
				node.referenced = False
			for spring in self.springs:
				if not spring.broken:
					value = getattr(spring, break_variable)
					value = value[0] if isinstance(value, list) else value
					spring.broken = SpringNetwork.relops[relop](value, threshold)
					if not spring.broken:
						spring.node_start.referenced = True
						spring.node_end.referenced = True
					else:
						any_new_breaks = True
			if any_new_breaks:
				self.connect_springs_nodes()


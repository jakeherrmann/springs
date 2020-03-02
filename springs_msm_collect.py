import sys
import random
import math
import statistics
import shutil
import struct
from pathlib import Path
import numpy as np
import springs_interface as spr 
import springs_biology as bio

def main(argv):

	#
	show_displays = False
	save_displays = False
	job_name = ''
	batch_name = ''
	if len(argv)>=1:
		job_name = argv[0]
	if len(argv)>=2:
		batch_name = argv[1]

	# file_format_nodes_new = spr.Node.get_file_format(2, 'double', use_solver_format=False)
	# file_format_nodes_old = spr.Node.get_file_format(2, 'double', use_solver_format=False)
	# for var in file_format_nodes_old.variables:
		# if var.name == 'fixed':
			# var.num = 1
			# var.py_type = bool
			
	# file_format_boundaries_new = spr.Boundary.get_file_format(2)
	# file_format_boundaries_old = spr.Boundary.get_file_format(2)
	# for var in file_format_boundaries_old.variables:
		# if var.name == 'fixed':
			# var.num = 1
			# var.py_type = bool
	
	# num_nodes = 5586
	# num_boundaries = 4
	# nodes = [ spr.Node(2) for count in range(num_nodes) ]
	# boundaries = [ spr.Boundary(2) for count in range(num_boundaries) ]
	# for path_job in sorted((Path('.')/'..'/'SAVE'/batch_name).glob('msm*_{:04d}'.format(int(job_name)))):
		# for path_stretch in sorted(path_job.glob( 'STRETCH_*' )):
			# #
			# file_name_old = path_stretch/'network_nodes.dat'
			# file_name_new = path_stretch/'network_nodes.new'
			# print( file_name_old )
			# file_format_nodes_old.read_binary_file(file_name_old, nodes)
			# for node in nodes:
				# node.fixed = [False]*2
			# file_format_nodes_new.write_binary_file(file_name_new, nodes)
			# #
			# file_name_old = path_stretch/'network_boundaries.dat'
			# file_name_new = path_stretch/'network_boundaries.new'
			# print( file_name_old )
			# file_format_boundaries_old.read_binary_file(file_name_old, boundaries)
			# for boundary in boundaries:
				# boundary.fixed = [False]*2
			# file_format_boundaries_new.write_binary_file(file_name_new, boundaries)
	# return
	
	num_nodes = 5586
	num_springs = 8273
	num_dimensions = 2
	
	format_str_nodes = '={:d}d'.format(num_nodes*num_dimensions)
	format_str_springs = '={:d}d'.format(num_springs)
	lung = bio.Lung()
	for path_job in sorted((Path('.')/'..'/'SAVE'/batch_name).glob('*_{:04d}'.format(int(job_name)))):
		file_name = str(path_job).replace('SAVE', 'COLLECT')+'.dat'
		Path(file_name).parent.mkdir(parents=True, exist_ok=True)
		with open(file_name,'wb') as file:	
			for path_stretch in sorted(path_job.glob( 'STRETCH_*' )):
				lung.net.read_spring_network(
					path_stretch,
					reinitialize=True)
				node_positions = [ n.position for n in lung.net.nodes ]
				node_positions = [ npi for np in node_positions for npi in np ]
				spring_stiffnesses = [ s.effective_stiffness() for s in lung.net.springs ]
				file.write( struct.pack(format_str_nodes, *node_positions) )
				file.write( struct.pack(format_str_springs, *spring_stiffnesses) )

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

	return

if __name__ == '__main__':
	main(sys.argv[1:])

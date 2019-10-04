import sys
import random
from pathlib import Path
import springs_interface as spr 
import springs_biology as bio

def main(argv):
	
	# create geometry for spring network and anatomical structures
	setup_type = 'hexagon_2D'

	if setup_type=='hexagon_2D':
		net, walls = spr.make_geom_hexagon_2D([8,8])
		net.precision = 'float'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		net.boundaries[0].fixed = True
		net.boundaries[1].fixed = True
		net.boundaries[2].fixed = True
		net.boundaries[3].fixed = True
		# net.boundaries[2].force_magnitudes = [+0.05] * len(net.boundaries[2].nodes)
		# net.boundaries[3].force_magnitudes = [+0.05] * len(net.boundaries[3].nodes)
		# for spring in net.springs:
		# 	spring.rest_length *= max(0.1, random.gauss(1.0,0.4))
		for spring in net.springs:
			# # TEST_01
			# spring.rest_length *= 0.7
			# k_mod = max(0.1, random.gauss(2.0,1.0))
			# # TEST_02
			# spring.rest_length *= 0.7
			# k_mod = max(0.1, random.gauss(2.0,0.2))
			# # TEST_03
			# spring.rest_length *= 0.147
			# k_mod = max(0.1, random.gauss(2.0,0.2))
			# # TEST_04
			# spring.rest_length *= 0.9
			# k_mod = max(0.1, random.gauss(2.0,0.2))
			# # TEST_05
			# spring.rest_length *= 1.1
			# k_mod = max(0.1, random.gauss(2.0,0.2))
			# # TEST_06
			# spring.rest_length *= 1.1
			# k_mod = max(0.1, random.gauss(2.0,1.0))
			# TEST_07
			spring.rest_length *= 1.0
			k_mod = max(0.1, random.gauss(2.0,1.0))
			# # TEST_08
			# spring.rest_length *= 0.95
			# k_mod = max(0.1, random.gauss(2.0,1.0))
			# # TEST_09
			# spring.rest_length *= 1.0
			# k_mod = max(0.1, random.gauss(2.0,1.0))
			spring.stiffness_tension = [ k*k_mod for k in spring.stiffness_tension ]

	elif setup_type=='truncoct_3D':
		net, walls = spr.make_geom_truncoct_3D([2,2,2])
		net.precision = 'float'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		net.boundaries[0].fixed = True
		net.boundaries[1].fixed = True
		net.boundaries[2].force_magnitudes = [+0.05] * len(net.boundaries[2].nodes)
		net.boundaries[3].force_magnitudes = [+0.05] * len(net.boundaries[3].nodes)
		net.boundaries[4].force_magnitudes = [+0.05] * len(net.boundaries[4].nodes)
		net.boundaries[5].force_magnitudes = [+0.05] * len(net.boundaries[5].nodes)
		for spring in net.springs:
			spring.rest_length *= max(0.1, random.gauss(1.0,0.4))

	elif setup_type=='file':
		net = spr.SpringNetwork()
		net.read_spring_network(Path('.')/'..'/'TEST', reinitialize=True)
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		walls = None

	else:
		print( 'UKNOWN SETUP' )
		return

	#
	lung = bio.Lung(
		spring_network=net,
		walls=walls,
		spring_break_variable=None,
		spring_break_threshold=4.0)
	lung.add_fibroblast_every_spring()
	current_stretch = 1.0
	stretches = [1.0]
	stretches += bio.stretch_triangle(stretch_amplitude=1.2, num_stretches_halfway=50)
	time_cycle = 5.0
	time_step = time_cycle / len(stretches)
	num_cycles = 50
	lung.stretch(1.0, dimensions='all', boundary_indexes=None)
	spr.display(lung.net,
		color_variable='stiffness_tension',
		color_range=(-0.05,5.0),
		ax_lims=((0.0,36.0),(0.0,36.0)),
		delay=None,
		save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(0),
		show=False)

	iter_total = 0
	for iter_cycle in range(num_cycles):
		for iter_stretch, stretch in enumerate(stretches):
			iter_total += 1
			current_stretch *= stretch ;
			print( '{:07.2f}% stretch, {:06d} springs remaining'.format(100.0*current_stretch, [ s.broken for s in lung.net.springs ].count(False)) )
			lung.stretch(stretch, dimensions='all', boundary_indexes=None)
			lung.save( Path('.')/'..'/'STRETCH_{:04d}'.format(iter_total) )
			lung.load( Path('.')/'..'/'STRETCH_{:04d}'.format(iter_total) )
			lung.agent_actions( time_step )
		# spr.display(lung.net,
		# 	color_variable='strain',
		# 	color_range=(-0.05,1.0),
		# 	ax_lims=((0.0,36.0),(0.0,36.0)),
		# 	delay=0.1,
		# 	save_file_name=None,
		# 	show=True)
		spr.display(lung.net,
			color_variable='stiffness_tension',
			color_range=(-0.05,5.0),
			ax_lims=((0.0,36.0),(0.0,36.0)),
			delay=None,
			save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(iter_cycle+1),
			show=False)

if __name__ == '__main__':
	main(sys.argv[1:])

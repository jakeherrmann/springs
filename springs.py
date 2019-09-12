import sys
import random
from pathlib import Path
import springs_interface as spr 
import springs_biology as bio

def main(argv):
	
	# create geometry for spring network and anatomical structures
	setup_type = 'hexagon_2D'

	if setup_type=='hexagon_2D':
		net, walls = spr.make_geom_hexagon_2D([18,18])
		net.precision = 'float'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		net.boundaries[0].fixed = True
		net.boundaries[1].fixed = True
		net.boundaries[2].fixed = True
		net.boundaries[3].fixed = True
		# net.boundaries[2].force_magnitudes = [+0.05] * len(net.boundaries[2].nodes)
		# net.boundaries[3].force_magnitudes = [+0.05] * len(net.boundaries[3].nodes)
		for spring in net.springs:
			spring.rest_length *= max(0.1, random.gauss(1.0,0.4))

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
		spring_break_variable='strain',
		spring_break_threshold=1.0)
	current_stretch = 1.0
	stretches = [1.0]
	stretches += bio.stretch_triangle(stretch_amplitude=2.0, num_stretches_halfway=10)
	for count, stretch in enumerate(stretches):
		current_stretch *= stretch ;
		print( '{:07.2f}% stretch, {:06d} springs remaining'.format(100.0*current_stretch, [ s.broken for s in lung.net.springs ].count(False)) )
		lung.stretch(stretch, dimensions='all', boundary_indexes=None)
		lung.save( Path('.')/'..'/'STRETCH_{:04d}'.format(count) )
		lung.load( Path('.')/'..'/'STRETCH_{:04d}'.format(count) )
		spr.display(lung.net,
			color_variable='strain',
			color_range=(-0.05,1.0),
			ax_lims=((0.0,80.0),(0.0,80.0)),
			delay=0.1,
			save_file_name=None,
			show=True)

if __name__ == '__main__':
	main(sys.argv[1:])

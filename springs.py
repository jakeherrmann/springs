import sys
import random
from pathlib import Path
import springs_interface as spr 

def main(argv):
	#
	setup_type = 'hexagon_2D'
	#
	if setup_type == 'hexagon_2D':
		net = spr.make_geom_hexagon_2D([16,16])
		net.precision = 'float'
	elif setup_type == 'file':
		net = spr.SpringNetwork()
		net.read_spring_network(Path('.')/'..'/'TEST')
	else:
		print( 'UKNOWN SETUP' )
		return
	#
	net.boundaries[0].fixed = True
	net.boundaries[1].fixed = True
	net.boundaries[2].force_magnitudes = [+0.2] * len(net.boundaries[2].nodes)
	net.boundaries[3].force_magnitudes = [+0.2] * len(net.boundaries[3].nodes)
	#
	for spring in net.springs:
		spring.rest_length = max(0.1, random.gauss(1.0,0.4))
	#
	net.apply_boundary_conditions()
	total_stretch = 3.0
	num_stretches = 100
	delta = (total_stretch-1.0)/num_stretches
	stretch_increments = [1.0] + [ (1.0+(n+1)*delta)/(1.0+n*delta) for n in range(num_stretches) ]
	current_stretch = 1.0
	for stretch in stretch_increments:
		current_stretch *= stretch ;
		print( '{:07.2f}% stretch, {:06d} springs remaining'.format(100.0*current_stretch, len(net.springs)) )
		net.dir_input  = Path('.')/'..'/'TEST_RESULT'/'INPUT'
		net.dir_output  = Path('.')/'..'/'TEST_RESULT'/'OUTPUT'
		# net.dir_output = Path('.')/'..'/'TEST_RESULT'/'STRAIN_{:06.0f}'.format(10000.0*(current_stretch-1.0))
		net.apply_stretch(stretch, 0)
		net.solve()
		net.calc_spring_force()
		img_file_name = (Path('.')/'..'/'TEST_RESULT'/'STRAIN_{:06.0f}'.format(10000.0*(current_stretch-1.0))).with_suffix('.png')
		spr.display(net,
			color_variable='force',
			color_range=(0.0,1.0),
			ax_lims=((0.0,80.0), (-15.0,+45.0)),
			delay=0.1,
			save_file_name=img_file_name,
			show=False)
		net.break_spring_strain(1.0)

if __name__ == '__main__':
	main(sys.argv[1:])

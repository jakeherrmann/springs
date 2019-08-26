import sys
from pathlib import Path
import springs_interface as spr 

def main(argv):
	net = spr.SpringNetwork()
	net.read_spring_network( Path('.') / '..' / 'TEST' )

	net.dir_input  = Path('.') / '..' / 'TEST_INPUT'
	net.dir_output = Path('.') / '..' / 'TEST_OUTPUT'
	net.solve()
	strains = [0.01] * 180
	total_stretch = 1.0 ;
	print( 'NUMBER OF SPRINGS REMAINING:' )
	print( net.num_springs )
	for strain in strains:
		net.stretch(strain, 0)
		total_stretch *= 1.0+strain ;
		net.dir_input  = Path('.') / '..' / 'TEST_{:06.0f}_INPUT'.format(10000.0*(total_stretch-1.0))
		net.dir_output = Path('.') / '..' / 'TEST_{:06.0f}_OUTPUT'.format(10000.0*(total_stretch-1.0))
		net.solve()
		net.calc_spring_force()
		net.break_spring_force(0.3)
		print( net.num_springs )

if __name__ == '__main__':
	main(sys.argv[1:])

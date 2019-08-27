import sys
from pathlib import Path
import springs_interface as spr 

def main(argv):
	net = spr.SpringNetwork()
	net.read_spring_network( Path('.') / '..' / 'TEST' )
	strains = [0.0] + ([0.01]*120)
	total_stretch = 1.0
	for strain in strains:
		total_stretch *= 1.0+strain ;
		print( '{:07.2f}% stretch, {:06d} springs remaining'.format(100.0*total_stretch, len(net.springs)) )
		net.dir_input  = Path('.') / '..' / 'TEST_RESULT' / 'INPUT'
		net.dir_output = Path('.') / '..' / 'TEST_RESULT' / 'STRETCH_{:06.0f}'.format(10000.0*(total_stretch-1.0))
		net.stretch(strain, 0)
		net.stretch(strain, 1)
		net.stretch(strain, 2)
		net.solve()
		net.calc_spring_force()
		net.break_spring_strain(4.0)

if __name__ == '__main__':
	main(sys.argv[1:])

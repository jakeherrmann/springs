import sys
import os
from pathlib import Path

def main(argv):
	dir_solver = Path('.') / 'springs_solver'
	dir_solver_obj = dir_solver / 'obj'
	os.system( 'mkdir ' + str(dir_solver_obj) )
	os.system( 'make -C ' + str(dir_solver) + ' clean' )
	os.system( 'make -C ' + str(dir_solver) )

if __name__ == '__main__':
	main(sys.argv[1:])

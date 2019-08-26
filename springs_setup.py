import sys
import os
from pathlib import Path

def main(argv):
	dir_springs_solver = Path('.') / 'springs_solver'
	os.system('make -C ./springs_solver clean')
	os.system('make -C ./springs_solver')

if __name__ == '__main__':
	main(sys.argv[1:])

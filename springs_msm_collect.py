import sys
import random
import math
import statistics
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
	if len(argv)>=2:
		job_name = argv[1]
	if len(argv)>=3:
		batch_name = argv[2]

	print( argv )
	print( list( (Path('.')/'SAVE'/batch_name).glob('*_{:04d}'.format(int(job_name))) ) )

	return

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	lung = bio.Lung()
	lung.net.read_spring_network(
		'',
		reinitialize=True)
	time_cycle = 5.0
	num_forces = 2
	num_cycles = 200
	iter_total = 0
	save_folder_name = 'msm_{:d}breath_{:d}D_force{:04.0f}-{:04.0f}'.format(
		num_cycles,
		lung.net.num_dimensions,
		1000*force_min,
		1000*force_max)
	save_folder_name += '_{:04d}'.format(int(job_name))
	lung.save( Path('.')/'..'/'SAVE'/batch_name/save_folder_name/'STRETCH_{:04d}'.format(iter_total) )

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

	return

if __name__ == '__main__':
	main(sys.argv[1:])

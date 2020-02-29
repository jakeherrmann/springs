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
	if len(argv)>=1:
		job_name = argv[0]
	if len(argv)>=2:
		batch_name = argv[1]

	lung = bio.Lung()
	for path_job in (Path('.')/'..'/'SAVE'/batch_name).glob('*_{:04d}'.format(int(job_name))):
		node_positions = []
		spring_stiffnesses = []
		for path_stretch in path_job.glob( 'STRETCH_*' ):
			lung.net.read_spring_network(
				path_stretch,
				reinitialize=True)
			node_positions.append( [ n.position for n in lung.net.nodes ] )
			spring_stiffnesses.append( [ s.effective_stiffness() for s in lung.net.springs ] )
			print( node_positions , spring_stiffnesses )

	return

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

	return

if __name__ == '__main__':
	main(sys.argv[1:])

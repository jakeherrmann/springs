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
	force_min   = 0.000
	force_delta = 0.025
	if len(argv)>=2:
		force_min   = float(argv[0])
		force_delta = float(argv[1])
	force_max = force_min + force_delta

	# create geometry for spring network and anatomical structures
	setup_type = 'hexagon_2D'

	if setup_type=='hexagon_2D':
		net = spr.make_geom_hexagon_2D([28,24]) ; fix_node = None #282 is top left
		# net = spr.make_geom_hexagon_2D([56,48])
		net.precision = 'double' #'float'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		# net.boundaries[0].fixed = True
		# net.boundaries[1].fixed = True
		# net.boundaries[2].fixed = True
		# net.boundaries[3].fixed = True
		net.boundaries[0].force_magnitudes = [+0.03] * len(net.boundaries[0].nodes)
		net.boundaries[1].force_magnitudes = [+0.03] * len(net.boundaries[1].nodes)
		net.boundaries[2].force_magnitudes = [+0.03] * len(net.boundaries[2].nodes)
		net.boundaries[3].force_magnitudes = [+0.03] * len(net.boundaries[3].nodes)
		# add some heterogeneity
		for spring in net.springs:
			spring.rest_length *= max(0.1, random.gauss(1.0,0.2))
			k_mod = max(0.1, random.gauss(1.0,0.2))
			spring.stiffness_tension = [ k*k_mod for k in spring.stiffness_tension ]

	elif setup_type=='truncoct_3D':
		net = spr.make_geom_truncoct_3D([6,6,6], split_walls_into_triangles=False)
		net.precision = 'double'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		# net.boundaries[0].fixed = True
		# net.boundaries[1].fixed = True
		net.boundaries[0].force_magnitudes = [+0.03] * len(net.boundaries[0].nodes)
		net.boundaries[1].force_magnitudes = [+0.03] * len(net.boundaries[1].nodes)
		net.boundaries[2].force_magnitudes = [+0.03] * len(net.boundaries[2].nodes)
		net.boundaries[3].force_magnitudes = [+0.03] * len(net.boundaries[3].nodes)
		net.boundaries[4].force_magnitudes = [+0.03] * len(net.boundaries[4].nodes)
		net.boundaries[5].force_magnitudes = [+0.03] * len(net.boundaries[5].nodes)
		# add some heterogeneity
		for spring in net.springs:
			spring.rest_length *= max(0.1, random.gauss(1.0,0.2))
			k_mod = max(0.1, random.gauss(1.0,0.2))
			spring.stiffness_tension = [ k*k_mod for k in spring.stiffness_tension ]

	elif setup_type=='file':
		net = spr.SpringNetwork()
		net.read_spring_network(Path('.')/'..'/'TEST', reinitialize=True)
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'

	else:
		print( 'UKNOWN SETUP' )
		return

	if not any( node.fixed for node in net.nodes ):
		if fix_node is not None:
			net.nodes[fix_node].fixed = True

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	net.solver_algorithm = 'steepest'
	net.solver_use_sum_net_force = True
	net.solver_tolerance_change_energy = 0.0
	net.solver_tolerance_sum_net_force = 1.0E-24
	net.solver_num_iter_print = 1000

	# net.solver_algorithm = 'newton'
	# net.solver_use_sum_net_force = True
	# net.solver_use_numerical_hessian = False
	# net.solver_tolerance_change_energy = 0.0
	# net.solver_tolerance_sum_net_force = 1.0E-24
	# net.solver_num_iter_print = 10

	# net.solver_algorithm = 'anneal'
	# net.solver_tolerance_change_energy = 1.0E-24
	# net.solver_tolerance_sum_net_force = 1.0E-02
	# net.solver_num_iter_print = 2000

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	using_stretch_profile = False
	using_forcing_profile = True

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	if using_forcing_profile:
		lung = bio.Lung(
			spring_network=net,
			spring_break_variable=None,
			spring_break_threshold=4.0)
		x = 70.0
		ax_lims = (
			(x*(1.0-1.8)/2.0,x*(1.0+1.8)/2.0),
			(x*(1.0-1.8)/2.0,x*(1.0+1.8)/2.0),
			(x*(1.0-1.8)/2.0,x*(1.0+1.8)/2.0))

		# note, this adds a fibroblast to every WALL, i.e., both sides of the same septal wall, i.e., two fibroblasts for each structure.
		lung.add_fibroblast_every_wall()

		# # for each wall, add opposite side of the septal wall to adjacent_walls
		# # this allows agents to move through walls
		# for wall in lung.walls:
		# 	for other_wall in lung.walls:
		# 		if wall is not other_wall:
		# 			if wall.structure is other_wall.structure:
		# 				if wall not in other_wall.adjacent_walls:
		# 					other_wall.adjacent_walls.append( wall )
		
		time_cycle = 5.0
		num_forces = 2
		num_cycles = 200
		iter_total = 0
		save_folder_name = 'msm_{:d}breath_{:d}D_hetero_small_force{:03.0f}-{:03.0f}'.format(
			num_cycles,
			lung.net.num_dimensions,
			1000*force_min,
			1000*force_max)
		if ( Path('.')/'..'/save_folder_name ).exists():
			print('Folder for this simulation already exists.  Exiting.')
			print(' ')
		else:
			for b in net.boundaries:
				b.force_magnitudes = [force_min] * len(b.nodes)
			lung.stretch(1.0, dimensions=None, boundary_indexes=None)
			if show_displays or save_displays:
				bio.display(lung,
					color_variable='stiffness_tension',
					color_range=(-0.05,20.0),
					ax_lims=ax_lims,
					delay=None,
					save_file_name=Path('.')/'..'/'breath_{:03d}.png'.format(0) if save_displays else None,
					show=show_displays)
			lung.save( Path('.')/'..'/save_folder_name/'STRETCH_{:04d}'.format(iter_total) )
			for iter_cycle in range(num_cycles):
				print(' ')
				print('{:03d}'.format(iter_cycle))
				# inspiration
				iter_total += 1
				for b in net.boundaries:
					b.force_magnitudes = [force_max] * len(b.nodes)
				lung.stretch(1.0, dimensions=None, boundary_indexes=None)
				lung.save( Path('.')/'..'/save_folder_name/'STRETCH_{:04d}'.format(iter_total) )
				spring_strains_in = [ spring.strain for spring in lung.net.springs ]
				current_stretch = sum(spring_strains_in) / len(spring_strains_in)
				print( ('IN:  '
						'{:07.2f}% average strain, '
						'{:06d} springs remaining, ').format(
						100.0*current_stretch,
						[ s.broken for s in lung.net.springs ].count(False)))
				# expiration
				iter_total += 1
				for b in net.boundaries:
					b.force_magnitudes = [force_min] * len(b.nodes)
				lung.stretch(1.0, dimensions=None, boundary_indexes=None)
				lung.save( Path('.')/'..'/save_folder_name/'STRETCH_{:04d}'.format(iter_total) )
				spring_strains_ex = [ spring.strain for spring in lung.net.springs ]
				current_stretch = sum(spring_strains_ex) / len(spring_strains_ex)
				print( ('EX:  '
						'{:07.2f}% average strain, '
						'{:06d} springs remaining, ').format(
						100.0*current_stretch,
						[ s.broken for s in lung.net.springs ].count(False)))
				# compute average strain and strain energy rate
				spring_strains_mean  = [ 0.5*(i+e) for i, e in zip(spring_strains_in, spring_strains_ex) ]
				spring_strains_delta = [  abs(i-e) for i, e in zip(spring_strains_in, spring_strains_ex) ]
				spring_strain_energy_rates = [ ((i**2)-(e**2))/time_cycle for i, e in zip(spring_strains_in, spring_strains_ex) ]
				spring_strain_energy_rates = [ r*(s.rest_length**2)*s.stiffness_tension[0] for s, r in zip(lung.net.springs, spring_strain_energy_rates) ]
				print( ('{:07.3f} mean average strain, '
						'{:07.5f} mean average strain energy rate').format(
						sum(spring_strains_mean)/len(spring_strains_mean),
						sum(spring_strain_energy_rates)/len(spring_strain_energy_rates)))

				# halt if solver fails
				if any( [ math.isnan(spring.strain) for spring in lung.net.springs ] ):
					print( 'FOUND NAN IN SOLUTION' )
					return

				# assign inputs to fibroblasts
				for agent in lung.agents:
					agent.strain = statistics.mean( [ spring_strains_mean[si] for si in agent.wall.structure.springs_indexes ] )
					agent.strain_energy_rate = statistics.mean( [ spring_strain_energy_rates[si] for si in agent.wall.structure.springs_indexes ] )
				lung.agent_actions( time_cycle )
				lung.net.break_spring('stiffness_tension', 0.0, relop='<=')
				if show_displays or save_displays:
					if (iter_cycle+1) % 5 is 0:
						bio.display(lung,
							color_variable='stiffness_tension',
							color_range=(-0.05,20.0),
							ax_lims=ax_lims,
							delay=None,
							save_file_name=Path('.')/'..'/'breath_{:03d}.png'.format(iter_cycle+1) if save_displays else None,
							show=show_displays)

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

	return

if __name__ == '__main__':
	main(sys.argv[1:])

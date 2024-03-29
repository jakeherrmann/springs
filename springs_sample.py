import sys
import random
import math
import statistics
from pathlib import Path
import springs_py_interface as spr 
import springs_py_biology as bio

def main(argv):
	
	# create geometry for spring network and anatomical structures
	setup_type = 'hexagon_2D'

	if setup_type=='hexagon_2D':
		net = spr.make_geom_hexagon_2D([28,24])
		net.precision = 'float'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		# net.boundaries[0].fixed = [True] * 2
		# net.boundaries[1].fixed = [True] * 2
		# net.boundaries[2].fixed = [True] * 2
		# net.boundaries[3].fixed = [True] * 2
		net.boundaries[0].force_magnitudes = [+0.03] * len(net.boundaries[0].nodes)
		net.boundaries[1].force_magnitudes = [+0.03] * len(net.boundaries[1].nodes)
		net.boundaries[2].force_magnitudes = [+0.03] * len(net.boundaries[2].nodes)
		net.boundaries[3].force_magnitudes = [+0.03] * len(net.boundaries[3].nodes)
		# add some heterogeneity
		for spring in net.springs:
			spring.rest_length *= max(0.1, random.gauss(1.0,0.2))
			k_mod = max(0.1, random.gauss(1.0,0.2))
			spring.force_length_parameters_tension = [ k*k_mod for k in spring.force_length_parameters_tension ]

	elif setup_type=='truncoct_3D':
		net = spr.make_geom_truncoct_3D([2,2,2], split_walls_into_triangles=True)
		net.precision = 'float'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		net.boundaries[0].fixed = [True] * 3
		net.boundaries[1].fixed = [True] * 3
		net.boundaries[2].force_magnitudes = [+0.03] * len(net.boundaries[2].nodes)
		net.boundaries[3].force_magnitudes = [+0.03] * len(net.boundaries[3].nodes)
		net.boundaries[4].force_magnitudes = [+0.03] * len(net.boundaries[4].nodes)
		net.boundaries[5].force_magnitudes = [+0.03] * len(net.boundaries[5].nodes)
		# add some heterogeneity
		for spring in net.springs:
			spring.rest_length *= max(0.1, random.gauss(1.0,0.2))
			k_mod = max(0.1, random.gauss(1.0,0.2))
			spring.force_length_parameters_tension = [ k*k_mod for k in spring.force_length_parameters_tension ]

	elif setup_type=='file':
		net = spr.SpringNetwork()
		net.read_spring_network(Path('.')/'..'/'TEST', reinitialize=True)
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'

	else:
		print( 'UKNOWN SETUP' )
		return

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	net.solver_algorithm = 'steepest'
	net.solver_objective = 'sumforce'
	net.solver_tolerance_change_objective = 1.0E-16
	net.solver_tolerance_sum_net_force = 1.0E-16
	net.solver_num_threads = 6
	net.solver_verbose = 1
	net.solver_parallel = 0

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	using_stretch_profile = False
	using_forcing_profile = True

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	if using_forcing_profile:
		lung = bio.Lung(
			spring_network=net,
			spring_break_variable=None,
			spring_break_threshold=4.0)
		# ax_lims = (
		# 	(42.5*(1.0-1.8)/2.0,42.5*(1.0+1.8)/2.0),
		# 	(42.5*(1.0-1.8)/2.0,42.5*(1.0+1.8)/2.0))
		ax_lims = (
			(42.5*(1.0-2.5)/2.0,42.5*(1.0+2.5)/2.0),
			(42.5*(1.0-2.5)/2.0,42.5*(1.0+2.5)/2.0))
		# ax_lims = (
		# 	(-4.0,12.0),
		# 	(-4.0,12.0),
		# 	(-4.0,12.0))

		# note, this adds a fibroblast to every WALL, i.e., both sides of the same septal wall, i.e., two fibroblasts for each structure.
		lung.add_fibroblast_every_wall()

		# for each wall, add opposite side of the septal wall to adjacent_walls
		# this allows agents to move through walls
		for wall in lung.walls:
			for other_wall in lung.walls:
				if wall is not other_wall:
					if wall.structure is other_wall.structure:
						if wall not in other_wall.adjacent_walls:
							other_wall.adjacent_walls.append( wall )
		
		# force_min = 0.425 ; force_max = 0.475
		force_min = 0.425 ; force_max = 0.525
		time_cycle = 5.0
		num_forces = 2
		num_cycles = 1 #600
		iter_total = 0
		for b in net.boundaries:
			b.force_magnitudes = [force_min] * len(b.nodes)
		lung.stretch(1.0, dimensions=None, boundary_indexes=None)
		bio.display(lung,
			color_variable=None,
			color_range=(-0.05,20.0),
			show_agents=False,
			ax_lims=ax_lims,
			delay=None,
			save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(0),
			dpi=300,
			show=False)
		lung.save( Path('.')/'..'/'SAVE'/'STRETCH_{:04d}'.format(iter_total) )
		for iter_cycle in range(num_cycles):
			print(' ')
			print('{:03d}'.format(iter_cycle))
			# inspiration
			iter_total += 1
			for b in net.boundaries:
				b.force_magnitudes = [force_max] * len(b.nodes)
			lung.stretch(1.0, dimensions=None, boundary_indexes=None)
			if (iter_cycle+1) % 10 == 0:
				lung.save( Path('.')/'..'/'STRETCH_{:04d}'.format(iter_total) )
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
			if (iter_cycle+1) % 10 == 0:
				lung.save( Path('.')/'..'/'STRETCH_{:04d}'.format(iter_total) )
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
			spring_strain_energy_rates = [ r*(s.rest_length**2)*s.force_length_parameters_tension[0] for s, r in zip(lung.net.springs, spring_strain_energy_rates) ]
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
			lung.net.break_spring('force_length_parameters_tension', 0.0, relop='<=')
			if (iter_cycle+1) % 10 == 0:
				bio.display(lung,
					color_variable='effective_spring_constant',
					color_range=(-0.05,5.0),
					show_agents=False,
					ax_lims=ax_lims,
					delay=None,
					save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(iter_cycle+1),
					dpi=300,
					show=False)
		bio.display(lung,
			color_variable='effective_spring_constant',
			color_range=(-0.05,5.0),
			show_agents=False,
			ax_lims=ax_lims,
			delay=None,
			save_file_name=Path('.')/'..'/'final.png'.format(iter_cycle+1),
			dpi=300,
			show=False)
		return

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	if using_stretch_profile:
		lung = bio.Lung(
			spring_network=net,
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
			color_variable=None,
			color_range=(-0.05,5.0),
			ax_lims=((15.0-20.0,15.0+20.0),(15.0-20.0,15.0+20.0)),
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
				lung.save( Path('.')/'..'/'SAVE'/'STRETCH_{:04d}'.format(iter_total) )
				lung.agent_actions( time_step )
			spr.display(lung.net,
				color_variable=None,
				color_range=(-0.05,5.0),
				ax_lims=((0.0,36.0),(0.0,36.0)),
				delay=None,
				save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(iter_cycle+1),
				show=False)
		return

if __name__ == '__main__':
	main(sys.argv[1:])

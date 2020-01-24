import sys
import random
import math
from pathlib import Path
import springs_interface as spr 
import springs_biology as bio

def main(argv):
	
	# create geometry for spring network and anatomical structures
	setup_type = 'hexagon_2D'

	if setup_type=='hexagon_2D':
		net, walls = spr.make_geom_hexagon_2D([28,24])
		net.precision = 'float'
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
			k_mod = max(0.1, random.gauss(2.0,1.0))
			spring.stiffness_tension = [ k*k_mod for k in spring.stiffness_tension ]

	elif setup_type=='truncoct_3D':
		net, walls = spr.make_geom_truncoct_3D([2,2,2])
		net.precision = 'float'
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		net.boundaries[0].fixed = True
		net.boundaries[1].fixed = True
		net.boundaries[2].force_magnitudes = [+0.03] * len(net.boundaries[2].nodes)
		net.boundaries[3].force_magnitudes = [+0.03] * len(net.boundaries[3].nodes)
		net.boundaries[4].force_magnitudes = [+0.03] * len(net.boundaries[4].nodes)
		net.boundaries[5].force_magnitudes = [+0.03] * len(net.boundaries[5].nodes)

	elif setup_type=='file':
		net = spr.SpringNetwork()
		net.read_spring_network(Path('.')/'..'/'TEST', reinitialize=True)
		net.dir_solver_input   = Path('.')/'..'/'SOLVER_DATA'/'INPUT'
		net.dir_solver_output  = Path('.')/'..'/'SOLVER_DATA'/'OUTPUT'
		walls = None

	else:
		print( 'UKNOWN SETUP' )
		return

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	using_stretch_profile = False
	using_forcing_profile = True

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	if using_forcing_profile:
		lung = bio.Lung(
			spring_network=net,
			walls=walls,
			spring_break_variable=None,
			spring_break_threshold=4.0)
		ax_lims = (
			(42.5*(1.0-1.8)/2.0,42.5*(1.0+1.8)/2.0),
			(42.5*(1.0-1.8)/2.0,42.5*(1.0+1.8)/2.0))

		lung.add_fibroblast_every_spring()
		force_min = 0.425 ; force_max = 0.475
		time_cycle = 5.0
		num_forces = 2
		num_cycles = 50
		iter_total = 0
		for b in net.boundaries:
			b.force_magnitudes = [force_min] * len(b.nodes)
		lung.stretch(1.0, dimensions=None, boundary_indexes=None)
		spr.display(lung.net,
			color_variable='stiffness_tension',
			color_range=(-0.05,20.0),
			ax_lims=ax_lims,
			delay=None,
			save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(0),
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
			spring_strain_energy_rates = [ r*(s.rest_length**2)*s.stiffness_tension[0] for s, r in zip(lung.net.springs, spring_strain_energy_rates) ]
			print( ('{:07.3f} mean average strain, '
					'{:07.5f} mean average strain energy rate').format(
					sum(spring_strains_mean)/len(spring_strains_mean),
					sum(spring_strain_energy_rates)/len(spring_strain_energy_rates)))
			# assign inputs to fibroblasts
			for agent in lung.agents:
				agent.strain = spring_strains_mean[ agent.spring_index ]
				agent.strain_energy_rate = spring_strain_energy_rates[ agent.spring_index ]
			lung.agent_actions( time_cycle )
			if (iter_cycle+1) % 5 is 0:
				spr.display(lung.net,
					color_variable='stiffness_tension',
					color_range=(-0.05,20.0),
					ax_lims=ax_lims,
					delay=None,
					save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(iter_cycle+1),
					show=False)
		return

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
	if using_stretch_profile:
		lung = bio.Lung(
			spring_network=net,
			walls=walls,
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
			color_variable='stiffness_tension',
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
				color_variable='stiffness_tension',
				color_range=(-0.05,5.0),
				ax_lims=((0.0,36.0),(0.0,36.0)),
				delay=None,
				save_file_name=Path('.')/'..'/'breath_{:02d}.png'.format(iter_cycle+1),
				show=False)
		return

if __name__ == '__main__':
	main(sys.argv[1:])

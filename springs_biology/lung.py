from .agent_fibroblast import Agent_Fibroblast
from .wall import Wall
from .alveolus import Alveolus

import springs_interface as spr
import numpy as np
import scipy.spatial

class Lung:
	def __init__(self,
			spring_network=None,
			agents=None,
			spring_break_variable=None,
			spring_break_threshold=0.0):
		if spring_network is None:
			self.net = spr.SpringNetwork()
		else:
			self.net = spring_network
		self.agents = agents if agents is not None else []
		self.spring_break_variable  = spring_break_variable
		self.spring_break_threshold = spring_break_threshold
		#
		self.alveoli = []
		self.walls = []
		if self.net is not None:
			self.build_lung_structure()

	def save(self, save_dir):
		self.net.write_spring_network(save_dir)

	def load(self, load_dir, reinitialize=False):
		self.net.read_spring_network(load_dir, reinitialize=reinitialize)
		if reinitialize:
			self.build_lung_structure()

	def stretch(self, stretch_increment, dimensions='all', boundary_indexes='all'):
		self.net.apply_stretch(stretch_increment, dimensions=dimensions, boundary_indexes=boundary_indexes)
		self.net.solve()
		self.net.calc_spring_force()
		if self.spring_break_variable is not None and self.spring_break_threshold is not None:
			self.net.break_spring(self.spring_break_variable, self.spring_break_threshold)

	def agent_actions(self, time_step):
		for agent in self.agents:
			agent.do_actions(time_step)

	def add_fibroblast_every_wall(self):
		for wall in self.walls:
			self.agents.append( Agent_Fibroblast(parent=self, wall=wall) )

	def build_lung_structure(self):
		# identify alveoli and 2-sided walls from structures and structure_groups
		self.alveoli = [ Alveolus() for sg in self.net.structure_groups ]
		self.walls.clear()
		for structure_group, alveolus in zip(self.net.structure_groups, self.alveoli):
			for structure in structure_group.structures:
				self.walls.append( Wall(structure=structure, alveolus=alveolus) )
				alveolus.walls.append( self.walls[-1] )
		for wall_index, wall in enumerate(self.walls):
			wall.index = wall_index
		# construct adjacency of walls within each alveolus
		for wall in self.walls:
			wall.adjacent_walls.clear()
		for wall in self.walls:
			for other_wall in wall.alveolus.walls:
				if wall is not other_wall:
					if wall.structure is other_wall.structure:
						other_wall.adjacent_walls.append( wall )
					elif any( n in other_wall.structure.nodes for n in wall.structure.nodes ):
						other_wall.adjacent_walls.append( wall )




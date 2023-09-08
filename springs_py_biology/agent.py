

class Agent:
	def __init__(self, parent=None, wall=None):
		self.parent = parent
		self.wall = wall
		if self.wall is not None:
			self.wall.agents.append( self )
		self.actions = None

	def do_actions(self, time_step, rate_multiplier=1.0):
		if self.actions is not None:
			for action in self.actions:
				action(time_step, rate_multiplier)
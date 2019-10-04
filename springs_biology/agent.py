

class Agent:
	def __init__(self):
		self.actions = None

	def do_actions(self, time_step):
		if self.actions is not None:
			for action in self.actions:
				action(time_step)
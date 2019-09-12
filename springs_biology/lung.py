
import springs_interface as spr

class Lung:
	def __init__(self,
			spring_network=None,
			walls=None,
			spring_break_variable=None,
			spring_break_threshold=0.0):
		self.net = spring_network
		self.walls = walls
		self.spring_break_variable  = spring_break_variable
		self.spring_break_threshold = spring_break_threshold
		self.alveoli = None
		self.fibroblasts = None

	def save(self, save_dir):
		self.net.write_spring_network(save_dir)
		file_name = save_dir/'lung_walls.dat'
		file_format = spr.Structure.get_file_format()
		file_format.write_binary_file(file_name, self.walls)

	def load(self, load_dir):
		self.net.read_spring_network(load_dir)
		file_name = load_dir/'lung_walls.dat'
		file_format = spr.Structure.get_file_format()
		file_format.read_binary_file(file_name, self.walls)

	def stretch(self, stretch_increment, dimensions='all', boundary_indexes='all'):
		self.net.apply_stretch(stretch_increment, dimensions=dimensions, boundary_indexes=boundary_indexes)
		self.net.solve()
		self.net.calc_spring_force()
		if self.spring_break_variable is not None and self.spring_break_threshold is not None:
			self.net.break_spring(self.spring_break_variable, self.spring_break_threshold)
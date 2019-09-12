import struct

class FileVariable:
	def __init__(self, name, file_type, num, py_type, active=True):
		self.name = name
		self.file_type = file_type
		self.num = num
		self.py_type = py_type
		self.active = active

class FileFormat:
	def __init__(self):
		self.fixed_format = False
		self.variables = []

	def read_binary_file(self, file_name, obj_list):
		if self.fixed_format:
			format_str = '=' + ''.join([ v.file_type*v.num for v in self.variables if v.active ])
			num_bytes = struct.calcsize(format_str)
			with open(file_name,'rb') as file:
				file_iter = struct.iter_unpack(format_str, file.read())
				for obj, data_values in zip(obj_list, file_iter):
					ind = 0
					for v in self.variables:
						if v.active:
							if isinstance(v.py_type, list):
								setattr(obj, v.name, list(data_values[ind:(ind+v.num)]))
							else:
								setattr(obj, v.name, data_values[ind])
							ind += v.num

		else:
			with open(file_name,'rb') as file:
				for obj in obj_list:
					for v in self.variables:
						if v.active:
							if isinstance(v.py_type, list):
								if v.num<0:
									input_num = struct.unpack('=I', file.read(struct.calcsize('I')))[0]
									format_str = '=' + v.file_type*input_num
									data_values = struct.unpack(format_str, file.read(struct.calcsize(format_str)))
									setattr(obj, v.name, list(data_values))
								else:
									format_str = '=' + v.file_type*v.num
									data_values = struct.unpack(format_str, file.read(struct.calcsize(format_str)))
									setattr(obj, v.name, list(data_values))
							else:
								format_str = '=' + v.file_type*v.num
								data_values = struct.unpack(format_str, file.read(struct.calcsize(format_str)))[0]
								setattr(obj, v.name, data_values)

	def write_binary_file(self, file_name, obj_list, fixed_format=True):
		if self.fixed_format:
			format_str = '=' + ''.join([ v.file_type*v.num for v in self.variables if v.active ])
			with open(file_name,'wb') as file:
				for obj in obj_list:
					output_values = []
					for v in self.variables:
						if v.active:
							if isinstance(v.py_type, list):
								output_list = getattr(obj, v.name)
								while output_list and isinstance(output_list[0], list):
									output_list = [ xi for x in output_list for xi in x ]
								if output_list:
									output_values.extend( output_list )
							else:
								output_values.append( getattr(obj, v.name) )
					file.write( struct.pack(format_str, *output_values) )

		else:
			with open(file_name,'wb') as file:
				for obj in obj_list:
					output_values = []
					format_str = '='
					for v in self.variables:
						if v.active:
							if isinstance(v.py_type, list):
								output_list = getattr(obj, v.name)
								while output_list and isinstance(output_list[0], list):
									output_list = [ xi for x in output_list for xi in x ]
								if v.num<0:
									if output_list:
										output_num = len(output_list)
										output_values.append(output_num)
										output_values.extend(output_list)
										format_str += 'I' + v.file_type*output_num
									else:
										output_values.append(0)
										format_str += 'I'
								else:
									output_values.extend(output_list)
									format_str += v.file_type*v.num
							else:
								output_values.append( getattr(obj, v.name) )
								format_str += v.file_type*v.num
					file.write( struct.pack(format_str, *output_values) )

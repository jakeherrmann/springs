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
		self.variables = []

	def read_binary_file(self, file_name, obj_list):
		format_str = ''.join([ v.file_type*v.num for v in self.variables if v.active ])
		num_bytes = struct.calcsize(format_str)
		with open(file_name,'rb') as file:
			file_data = file.read()
			num_data = int( len(file_data) / num_bytes )
			file_iter = struct.iter_unpack(format_str,file_data)
			for data_index, data_values in enumerate(file_iter):
				byte_index = 0
				for v in self.variables:
					if isinstance(v.py_type, list):
						setattr(obj_list[data_index], v.name, list(data_values[ byte_index:(byte_index+v.num)]))
					else:
						setattr(obj_list[data_index], v.name, data_values[ byte_index ])
					byte_index += v.num

	def write_binary_file(self, file_name, obj_list):
		format_str = ''.join([ v.file_type*v.num for v in self.variables if v.active ])
		with open(file_name,'wb') as file:
			for obj in obj_list:
				output_values = []
				for v in self.variables:
					if isinstance(v.py_type, list):
						output_values += getattr(obj, v.name)
					else:
						output_values.append( getattr(obj, v.name) )
				file.write( struct.pack(format_str, *output_values) )
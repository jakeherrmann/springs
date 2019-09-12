import math

def stretch_ramp(stretch_amplitude, num_stretches):
	delta = (stretch_amplitude-1.0)/num_stretches
	stretch_increments = [ (1.0+(n+1)*delta)/(1.0+n*delta) for n in range(num_stretches) ]
	return stretch_increments

def stretch_triangle(stretch_amplitude, num_stretches_halfway):
	delta = (stretch_amplitude-1.0)/num_stretches_halfway
	stretch_increments = [ (1.0+(n+1)*delta)/(1.0+n*delta) for n in range(num_stretches_halfway) ]
	stretch_increments = stretch_increments[:] + [ 1.0/si for si in stretch_increments[::-1] ]
	return stretch_increments
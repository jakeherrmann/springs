
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mplcol

def display(spring_network, color_variable=None, color_range=None, ax_lims=None, delay=None, save_file_name=False, show=True):
	cmap = plt.get_cmap('inferno')
	if color_variable is None:
		c = [(0.0,0.0,0.0)] * len(spring_network.springs)
	else:
		if hasattr(spring_network.springs[0], color_variable):
			v = [ getattr(spring,color_variable) for spring in spring_network.springs ]
		else:
			v = [0.0] * len(spring_network.springs)
		if color_range is None:
			vmin = min(v)
			vmax = max(v)
			color_range = (vmin, vmax) if vmin!=vmax else (vmin,vmin+1.0)
		w = [ (vi-color_range[0])/(color_range[1]-color_range[0]) for vi in v ]
		w = [ max(0.0, min(1.0, wi)) for wi in w ]
		w = [ 1.0-wi for wi in w ]
		c = [ cmap(wi) for wi in w ]

	xyz_segments = [
		np.stack([
			np.asarray(spring.node_start_pointer.position),
			np.asarray(spring.node_end_pointer.position  )])
		for spring in spring_network.springs ]

	plt.ioff()
	if not plt.get_fignums():
		plt.gcf().set_size_inches(10.0, 6.0)
	ax = plt.gca()
	plt.cla()
	ax.set_aspect('equal')
	if ax_lims is not None:
		ax.set_xlim(ax_lims[0][0], ax_lims[0][1])
		ax.set_ylim(ax_lims[1][0], ax_lims[1][1])
	ax.add_collection( mplcol.LineCollection(xyz_segments, colors=c) )
	if save_file_name is not None:
		plt.savefig(save_file_name)
	if not show:
		plt.close(plt.gcf())
	else:
		if delay == float('inf'):
			plt.show()
		elif delay is not None:
			plt.draw()
			plt.pause(delay)

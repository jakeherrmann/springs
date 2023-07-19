
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mplcol

import scipy.spatial

import matplotlib.colors as mplc
import mpl_toolkits.mplot3d as mpl3d

def display(lung, show_slice=False, **kwargs):
	if lung.net.num_dimensions==2:
		display_2D(lung, **kwargs)
	elif lung.net.num_dimensions==3:
		if show_slice:
			display_3D_slice(lung, **kwargs)
		else:
			display_3D(lung, **kwargs)

def display_2D(lung, color_variable=None, color_range=None, ax_lims=None, delay=None, save_file_name=None, show_agents=True, show=True, dpi=300):
	spring_network = lung.net
	cmap = plt.get_cmap('inferno')
	if color_variable is None:
		colors = [(0.0,0.0,0.0)] * len(spring_network.springs)
	else:
		if hasattr(spring_network.springs[0], color_variable):
			if isinstance(getattr(spring_network.springs[0],color_variable), list):
				v = [ getattr(spring,color_variable)[0] for spring in spring_network.springs ]
			else:
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
		colors = [ cmap(wi) for wi in w ]

	segments = [
		np.stack([
			np.asarray(spring.node_start.position),
			np.asarray(spring.node_end.position  )])
		for spring in spring_network.springs ]

	segments = [ s for s, spring in zip(segments, spring_network.springs) if not spring.broken ]
	colors   = [ c for c, spring in zip(colors,   spring_network.springs) if not spring.broken ]

	plt.ioff()
	if not plt.get_fignums():
		plt.gcf().set_size_inches(10.0, 6.0)
	ax = plt.gca()
	plt.cla()
	ax.set_aspect('equal')
	if ax_lims is not None:
		ax.set_xlim(ax_lims[0][0], ax_lims[0][1])
		ax.set_ylim(ax_lims[1][0], ax_lims[1][1])
	else:
		s_max = np.asarray( [ s.max(0) for s in segments ] ).max(0)
		s_min = np.asarray( [ s.min(0) for s in segments ] ).min(0)
		ax.set_xlim(s_min[0], s_max[0])
		ax.set_ylim(s_min[1], s_max[1])
	ax.add_collection( mplcol.LineCollection(segments, colors=colors) )

	# alv_centroids = np.array([ np.mean( np.array( [ n.position for w in a.walls for n in w.structure.nodes ] ) ,axis=0) for a in lung.alveoli ])
	# plt.scatter(
	# 	alv_centroids[:,0],
	# 	alv_centroids[:,1],
	# 	marker='o',
	# 	edgecolor=(0.0,0.0,0.0),
	# 	facecolor='None')

	if show_agents and lung.agents:
		agent_positions = np.array( [ np.mean( np.array( [ n.position for n in a.wall.structure.nodes ] ) ,axis=0) for a in lung.agents ] )
		plt.scatter(
			agent_positions[:,0],
			agent_positions[:,1],
			marker='o',
			edgecolor=(0.0,0.0,0.0),
			facecolor=(0.0,0.0,0.0))
	
	if save_file_name is not None:
		plt.savefig(save_file_name, dpi=dpi)
	if not show:
		plt.close(plt.gcf())
	else:
		if delay == float('inf'):
			plt.show()
		elif delay is not None:
			plt.draw()
			plt.pause(delay)

def display_3D_slice(lung, color_variable=None, color_range=None, ax_lims=None, delay=None, save_file_name=None, show=True, structures=None, dpi=300):
	spring_network = lung.net
	plane_point  = np.array([2.0, 2.0, 5.0])
	plane_normal = np.array([0.0, 0.0, 1.0])
	plane_basis  = [
		np.array([1.0, 0.0, 0.0]),
		np.array([0.0, 1.0, 0.0])]
	cmap = plt.get_cmap('inferno')
	xyz_segments = []
	c = []
	for structure in structures:
		tri_points = np.stack([ np.asarray(n.position) for n in structure.nodes ])
		if tri_points.shape[0]==3:
			# plane-triangle intersection
			sideness = np.dot(tri_points-plane_point, plane_normal)
			above = np.where(sideness> 0.0)[0]
			below = np.where(sideness<=0.0)[0]
			if above.size==1:
				i0 = above[0]
				i1 = below[0]
				i2 = below[1]
			elif below.size==1:
				i0 = below[0]
				i1 = above[0]
				i2 = above[1]
			else:
				continue
			p0 = tri_points[i0,:]
			p1 = tri_points[i1,:]
			p2 = tri_points[i2,:]
			r01 = np.dot(plane_point-p0, plane_normal) / np.dot(p1-p0, plane_normal)
			r02 = np.dot(plane_point-p0, plane_normal) / np.dot(p2-p0, plane_normal)
			x01 = p0 + r01*(p1-p0)
			x02 = p0 + r02*(p2-p0)
			# convert to coordinates on 2D plane
			v01 = np.array([ np.dot(x01-plane_point, pb) for pb in plane_basis ])
			v02 = np.array([ np.dot(x02-plane_point, pb) for pb in plane_basis ])
			xyz_segments.append(np.stack([v01,v02]))
			# color
			if color_variable is None:
				ci = (0.0,0.0,0.0)
			else:
				if hasattr(structure.springs[0], color_variable):
					v = [ getattr(spring, color_variable) for spring in structure.springs ]
				else:
					v = [0.0] * len(structure.springs)
				if color_range is None:
					vmin = min(v)
					vmax = max(v)
					color_range = (vmin, vmax) if vmin!=vmax else (vmin,vmin+1.0)
				wi = sum(v)/len(v)
				wi = (wi-color_range[0])/(color_range[1]-color_range[0])
				wi = max(0.0, min(1.0, wi))
				wi = 1.0-wi
				ci = cmap(wi)
			c.append(ci)

	xyz_min = np.asarray( [ xyz.min(0) for xyz in xyz_segments ] ).min(0)
	xyz_segments = [ xyz-xyz_min for xyz in xyz_segments ]

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
		plt.savefig(save_file_name, dpi=dpi)
	if not show:
		plt.close(plt.gcf())
	else:
		if delay == float('inf'):
			plt.show()
		elif delay is not None:
			plt.draw()
			plt.pause(delay)

def display_3D(lung, color_variable=None, color_range=None, ax_lims=None, delay=None, save_file_name=False, show_agents=True, show=True, structures=None, dpi=300):
	spring_network = lung.net
	xyz_segments = [
		np.stack([
			np.asarray(spring.node_start.position),
			np.asarray(spring.node_end.position  )])
		for spring in spring_network.springs ]

	plt.ioff()
	if not plt.get_fignums():
		plt.gcf().set_size_inches(10.0, 6.0)
	ax = plt.gca()
	plt.cla()
	ax = plt.axes(
		projection='3d',
		proj_type='ortho')
	ax.set_position((0.0,0.0,1.0,1.0))
	#ax.set_aspect('equal')
	if ax_lims is not None:
		ax.set_xlim(ax_lims[0][0], ax_lims[0][1])
		ax.set_ylim(ax_lims[1][0], ax_lims[1][1])
		ax.set_zlim(ax_lims[2][0], ax_lims[2][1])
	else:
		xyz_max = np.asarray( [ xyz.max(0) for xyz in xyz_segments ] ).max(0)
		xyz_min = np.asarray( [ xyz.min(0) for xyz in xyz_segments ] ).min(0)
		ax.set_xlim(xyz_min[0], xyz_max[0])
		ax.set_ylim(xyz_min[1], xyz_max[1])
		ax.set_zlim(xyz_min[2], xyz_max[2])
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
	ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
	ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
	ax.xaxis._axinfo["grid"]['color'] = (1.0,1.0,1.0,0.0)
	ax.yaxis._axinfo["grid"]['color'] = (1.0,1.0,1.0,0.0)
	ax.zaxis._axinfo["grid"]['color'] = (1.0,1.0,1.0,0.0)
	if structures is not None:
		xyz_walls = [
			np.stack([ np.asarray(node.position) for node in structure.nodes ])
			for structure in structures ]
		cmap = plt.get_cmap('inferno')
		if color_variable is None:
			c = [(0.3,0.3,1.0)] * len(structures)
		else:
			if hasattr(spring_network.springs[0], color_variable):
				v = [ getattr(spring, color_variable) for spring in spring_network.springs ]
			else:
				v = [0.0] * len(spring_network.springs)
			if color_range is None:
				vmin = min(v)
				vmax = max(v)
				color_range = (vmin, vmax) if vmin!=vmax else (vmin,vmin+1.0)
			w = [ sum([ v[s] for s in w.springs ])/len(w.springs) for w in structures ]
			w = [ (wi-color_range[0])/(color_range[1]-color_range[0]) for wi in w ]
			w = [ max(0.0, min(1.0, wi)) for wi in w ]
			w = [ 1.0-wi for wi in w ]
			c = [ cmap(wi) for wi in w ]
		ax.add_collection3d(mpl3d.art3d.Poly3DCollection(
			xyz_walls,
			edgecolor=(0.0,0.0,0.0),
			linewidth=0.5,
			alpha=0.6,
			facecolor=c))
	else:
		cmap = plt.get_cmap('inferno')
		if color_variable is None:
			c = [(0.0,0.0,0.0)] * len(spring_network.springs)
		else:
			if hasattr(spring_network.springs[0], color_variable):
				if isinstance(getattr(spring_network.springs[0],color_variable), list):
					v = [ getattr(spring,color_variable)[0] for spring in spring_network.springs ]
				else:
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
		ax.add_collection3d(mpl3d.art3d.Line3DCollection(
			xyz_segments,
			color=c))

	# alv_centroids = np.array([ np.mean( np.array( [ n.position for w in a.walls for n in w.structure.nodes ] ) ,axis=0) for a in lung.alveoli ])
	# alv_centroids = np.array([ ac for ac in alv_centroids if not np.any(np.isnan(ac)) ])
	# ax.scatter(
	# 	alv_centroids[:,0],
	# 	alv_centroids[:,1],
	# 	alv_centroids[:,2],
	# 	marker='o',
	# 	edgecolor=(0.0,0.0,0.0),
	# 	facecolor='None')

	if show_agents and lung.agents:
		agent_positions = np.array( [ np.mean( np.array( [ n.position for n in a.wall.structure.nodes ] ) ,axis=0) for a in lung.agents ] )
		ax.scatter(
			agent_positions[:,0],
			agent_positions[:,1],
			agent_positions[:,2],
			marker='o',
			edgecolor=(0.0,0.0,0.0),
			facecolor=(0.0,0.0,0.0))

	if save_file_name is not None:
		plt.savefig(save_file_name, dpi=dpi)
	if not show:
		plt.close(plt.gcf())
	else:
		if delay == float('inf'):
			plt.show()
		elif delay is not None:
			plt.draw()
			plt.pause(delay)


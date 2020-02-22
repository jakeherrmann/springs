from ..node           import Node
from ..spring         import Spring
from ..boundary       import Boundary
from ..structure      import Structure
from ..structureGroup import StructureGroup
from ..springNetwork  import SpringNetwork
import math
import numpy as np
import scipy.spatial
# import matplotlib.pyplot as plt

def make_geom_hexagon_2D(geom_size=[1,1]):
	num_row = geom_size[1] + 3
	num_col = geom_size[0] + 4
	sqrt3 = math.sqrt(3.0)

	# obtain vertices and edges from a voronoi diagram
	# adjust points for unit hexagon side lengths
	xyz = np.indices((num_col,num_row), dtype=np.float64)
	xyz[0] *= 0.5*sqrt3
	xyz[1] += 0.5*np.mod( np.arange(0,num_col,1)+1 ,2)[:,np.newaxis]
	xyz *= sqrt3
	xyz = xyz.reshape((2,-1)).transpose()
	vor = scipy.spatial.Voronoi( xyz )
	v = vor.vertices
	e = vor.ridge_vertices
	r = vor.regions

	# scipy.spatial.voronoi_plot_2d(vor)
	# plt.show()
	
	# remove extraneous vertices
	small_number = 0.01
	e = [ ei for ei in e if all([ eij>=0 for eij in ei ]) ]
	e = [ ei for ei in e if all([ np.all(v[eij,:]>=0.0)                         for eij in ei ]) ]
	e = [ ei for ei in e if all([ v[eij,1]<=((num_row-1)*sqrt3)  + small_number for eij in ei ]) ]
	e = [ ei for ei in e if all([ v[eij,0]>0.5                   + small_number for eij in ei ]) ]
	e = [ ei for ei in e if all([ v[eij,0]<((num_col-1)*1.5-0.5) - small_number for eij in ei ]) ]
	e = [ ei for ei in e if any([ v[eij,1]>(0.5*sqrt3)           + small_number for eij in ei ]) ]
	e = [ ei for ei in e if any([ v[eij,1]<((num_row-1)*sqrt3)   - small_number for eij in ei ]) ]

	# separate boundary edges (to neighboring alveoli) from internal edges
	b = [
		[ ei for ei in e if np.any(np.isclose(v[ei,0],1.0                )) ],
		[ ei for ei in e if np.any(np.isclose(v[ei,0],(num_col-1)*1.5-1.0)) ],
		[ ei for ei in e if np.any(np.isclose(v[ei,1],0.5*sqrt3          )) ],
		[ ei for ei in e if np.any(np.isclose(v[ei,1],(num_row-1)*sqrt3  )) ]
		]
	for bi in b:
		e = [ ej for ej in e if ej not in bi ]

	# construct nodes
	nodes_ind = list(set([ eij for ei in e for eij in ei ]))
	nodes_position = v[ nodes_ind ,:]
	nodes_position -= nodes_position.min(axis=0)
	nodes = [ Node(num_dimensions=2) for count in range(len(nodes_ind)) ]
	for n in range(len(nodes)):
		nodes[n].position = list(nodes_position[n,:])

	# construct springs
	springs_nodes = [ [ nodes_ind.index(eij) for eij in ei ] for ei in e ]
	springs = [ Spring(node_start_index=sn[0], node_end_index=sn[1]) for sn in springs_nodes ]
	for spring in springs:
		spring.rest_length = np.linalg.norm( nodes_position[spring.node_end_index,:] - nodes_position[spring.node_start_index,:] )

	# construct wall structures
	structures = [ Structure(nodes_indexes=[s.node_start_index, s.node_end_index], springs_indexes=[i]) for i, s in enumerate(springs) ]

	# construct regions enclosed by wall structures
	regions_nodes_indexes = [ [ nodes_ind.index(ri) for ri in r if ri in nodes_ind ] for r in vor.regions ]
	nodes_regions = [ [] for n in nodes ]
	for region_index, region_nodes_indexes in enumerate(regions_nodes_indexes):
		for n in region_nodes_indexes:
			nodes_regions[n].append( region_index )
	regions_structures = [ [] for r in vor.regions ]
	for structure_index, structure in enumerate(structures):
		for node_index in structure.nodes_indexes:
			for region_index in nodes_regions[node_index]:
				if all( n in regions_nodes_indexes[region_index] for n in structure.nodes_indexes ):
					regions_structures[region_index].append( structure_index )
	structure_groups = [ StructureGroup(structures_indexes=rs) for rs in regions_structures if len(rs) > 1 ]

	# construct boundaries
	for bi in b:
		for bij in bi:
			if bij[1] in nodes_ind:
				bij.reverse()
	boundaries = [ Boundary(num_dimensions=2) for bi in b ]
	boundaries[0].outward_direction = [-1.0, 0.0]
	boundaries[1].outward_direction = [+1.0, 0.0]
	boundaries[2].outward_direction = [ 0.0,-1.0]
	boundaries[3].outward_direction = [ 0.0,+1.0]
	for boundary, bi in zip(boundaries, b):
		boundary.nodes_indexes = [ nodes_ind.index(bij[0]) for bij in bi ]
		start = [ bij[0] for bij in bi ]
		end   = [ bij[1] for bij in bi ]
		bv = v[end,:] - v[start,:]
		bv /= np.linalg.norm(bv, axis=1, keepdims=True)
		boundary.force_directions = bv.tolist()

	#
	net = SpringNetwork(num_dimensions=2)
	net.setup(
		nodes=nodes,
		springs=springs,
		structures=structures,
		structure_groups=structure_groups,
		boundaries=boundaries)

	return net

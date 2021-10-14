from ..node           import Node
from ..spring         import Spring
from ..boundary       import Boundary
from ..structure      import Structure
from ..structureGroup import StructureGroup
from ..springNetwork  import SpringNetwork
import math
import numpy as np
import scipy.spatial

def make_geom_truncoct_3D(geom_size=[5,3,4], split_walls_into_triangles=False):
	num_row = geom_size[1] + 3
	num_col = geom_size[0] + 3
	num_lay = geom_size[2] + 3

	# obtain vertices and edges from a voronoi diagram
	xyz = np.indices((2*num_col,num_row,num_lay), dtype=np.float64)
	xyz[0] *= 0.5
	xyz[1] += 0.5 * np.mod(np.arange(0,num_col*2,1)+1,2)[:,np.newaxis,np.newaxis]
	xyz[2] += 0.5 * np.mod(np.arange(0,num_col*2,1)+1,2)[:,np.newaxis,np.newaxis]
	xyz = xyz.reshape((3,-1)).transpose()
	vor = scipy.spatial.Voronoi( xyz )
	v = vor.vertices
	e = vor.ridge_vertices

	# remove extraneous vertices
	e = [ ei for ei in e if all([ eij>=0 for eij in ei ]) ]
	e = [ ei for ei in e if all([ np.all(v[eij,:]>=0.0       ) for eij in ei ]) ]
	e = [ ei for ei in e if all([ np.all(v[eij,:]>=(0.99*0.5)) for eij in ei ]) ]
	e = [ ei for ei in e if all([ v[eij,0]<=(num_col-1) for eij in ei ]) ]
	e = [ ei for ei in e if all([ v[eij,1]<=(num_row-1) for eij in ei ]) ]
	e = [ ei for ei in e if all([ v[eij,2]<=(num_lay-1) for eij in ei ]) ]

	# separate boundary edges (to neighboring alveoli) from internal edges
	b = [
		[ ei for ei in e if any([ v[eij,0]<1.0*0.99           for eij in ei ]) ],
		[ ei for ei in e if any([ v[eij,0]>(num_col-1.5)*1.01 for eij in ei ]) ],
		[ ei for ei in e if any([ v[eij,1]<1.0*0.99           for eij in ei ]) ],
		[ ei for ei in e if any([ v[eij,1]>(num_row-1.5)*1.01 for eij in ei ]) ],
		[ ei for ei in e if any([ v[eij,2]<1.0*0.99           for eij in ei ]) ],
		[ ei for ei in e if any([ v[eij,2]>(num_lay-1.5)*1.01 for eij in ei ]) ]
		]
	for bi in b:
		e = [ ej for ej in e if ej not in bi ]

	# adjust points for unit truncated octahedron side lengths
	v *= math.sqrt(8.0)

	# construct nodes
	nodes_ind = list(set([ eij for ei in e for eij in ei ]))
	nodes_position = v[ nodes_ind ,:]
	nodes_position -= nodes_position.min(axis=0)
	nodes = [ Node(num_dimensions=3) for count in range(len(nodes_ind)) ]
	for n in range(len(nodes)):
		nodes[n].position = nodes_position[n,:].tolist()

	# identify groups of nodes and springs that form walls
	walls_nodes = [ [ nodes_ind.index(eij) for eij in ei ] for ei in e ]
	springs_nodes = [ [wn[i-1], wn[i]] for wn in walls_nodes for i in range(len(wn)) ]
	walls_springs = [ range(len(wn)) for wn in walls_nodes ]
	for i in range(len(walls_springs)-1):
		count = walls_springs[i][-1] + 1
		walls_springs[i+1] = [ ws+count for ws in walls_springs[i+1] ]
	springs_nodes, unique_ind = np.unique(springs_nodes, axis=0, return_inverse=True)
	springs_nodes = springs_nodes.tolist()
	walls_springs = [ [ unique_ind[wsi] for wsi in ws ] for ws in walls_springs ]

	# identify groups of walls that form enclosed regions
	regions_walls = [ [] for r in vor.regions ]
	regions_nodes_indexes = [ [ nodes_ind.index(ri) for ri in r if ri in nodes_ind ] for r in vor.regions ]
	for wall_index, wall_nodes in enumerate(walls_nodes):
		for region_index, region_nodes_indexes in enumerate(regions_nodes_indexes):
			if all( n in region_nodes_indexes for n in wall_nodes ):
				regions_walls[region_index].append( wall_index )

	if split_walls_into_triangles:
		# for each wall, optionally add a central node and springs connecting to all other nodes
		tri_regions_walls = [ [] for r in vor.regions ]
		tri_walls_nodes = []
		tri_walls_springs = []
		for wi, (wn, ws) in enumerate(zip(walls_nodes, walls_springs)):
			centroid = np.mean( nodes_position[wn,:] ,axis=0)
			new_node = len(nodes)
			new_springs_nodes = [ [wni, new_node] for wni in wn ]
			new_springs = [ len(springs_nodes)+i for i in range(len(new_springs_nodes)) ]
			new_tri_nodes = [ [wn[i-1],wn[i],new_node] for i in range(len(wn)) ]
			new_tri_springs = [ [ws[i],new_springs[i-1],new_springs[i]] for i in range(len(ws)) ]
			new_tri_walls = range(len(tri_walls_nodes), len(tri_walls_nodes)+len(new_tri_nodes)-1)
			for region, tri_region in zip(regions_walls, tri_regions_walls):
				if wi in region:
					tri_region.extend( new_tri_walls )
			tri_walls_nodes.extend( new_tri_nodes )
			tri_walls_springs.extend( new_tri_springs )
			wn.append( new_node )
			ws.extend( new_springs )
			springs_nodes.extend( new_springs_nodes )
			nodes.append( Node(num_dimensions=3) )
			nodes[new_node].position = centroid.tolist()
		regions_walls = tri_regions_walls
		walls_nodes   = tri_walls_nodes
		walls_springs = tri_walls_springs

	# construct springs
	springs = [ Spring(node_start_index=sn[0], node_end_index=sn[1]) for sn in springs_nodes ]
	for spring in springs:
		p_end   = np.asarray(nodes[spring.node_end_index  ].position)
		p_start = np.asarray(nodes[spring.node_start_index].position)
		spring.rest_length = np.linalg.norm( p_end - p_start )

	# construct wall structures
	structures = [ Structure(nodes_indexes=wn, springs_indexes=ws) for wn, ws in zip(walls_nodes, walls_springs) ]

	# construct regions enclosed by wall structures
	structure_groups = [ StructureGroup(structures_indexes=rw) for rw in regions_walls if len(rw) > 1 ]

	# identify pairs of internal nodes and corresponding external boundary nodes
	b = [ [ bij for bij in bi if any([ bijk in nodes_ind for bijk in bij ]) ] for bi in b ]
	B = []
	for i, bi in enumerate(b):
		B.append([])
		for j, bij in enumerate(bi):
			for k, bijk in enumerate(bij):
				if bijk in nodes_ind:
					k_prev = k-1 if (k-1)>0        else len(bij)-1
					k_next = k+1 if (k+1)<len(bij) else 0
					if bij[k_prev] not in nodes_ind:
						B[i].append( (bijk, bij[k_prev]) )
					if k_next!=k_prev and bij[k_next] not in nodes_ind:
						B[i].append( (bijk, bij[k_next]) )

	# construct boundaries
	boundaries = [ Boundary(num_dimensions=3) for bi in b ]
	boundaries[0].outward_direction = [-1.0, 0.0, 0.0]
	boundaries[1].outward_direction = [+1.0, 0.0, 0.0]
	boundaries[2].outward_direction = [ 0.0,-1.0, 0.0]
	boundaries[3].outward_direction = [ 0.0,+1.0, 0.0]
	boundaries[4].outward_direction = [ 0.0, 0.0,-1.0]
	boundaries[5].outward_direction = [ 0.0, 0.0,+1.0]
	for boundary, Bi in zip(boundaries, B):
		B_nodes_ind = [ nodes_ind.index(Bij[0]) for Bij in Bi ] ;
		boundary.nodes_indexes = list(set(B_nodes_ind))
		start = [ Bij[0] for Bij in Bi ]
		end   = [ Bij[1] for Bij in Bi ]
		Bv = v[end,:] - v[start,:]
		Bv /= np.linalg.norm(Bv, axis=1, keepdims=True)
		nodes_force = np.zeros((len(nodes), 3), dtype=float)
		np.add.at( nodes_force , B_nodes_ind , Bv )
		boundary.force_directions = nodes_force[boundary.nodes_indexes,:].tolist()

	#
	net = SpringNetwork(num_dimensions=3)
	net.setup(
		nodes=nodes,
		springs=springs,
		structures=structures,
		structure_groups=structure_groups,
		boundaries=boundaries)

	return net

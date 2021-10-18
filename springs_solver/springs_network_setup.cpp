#include "springs.hpp"

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <map>

#include <metis.h>

#include <omp.h>
#include <stdlib.h>

#ifdef _WIN32
	#define FILESEP '\\'
#else
	#define FILESEP '/'
#endif

#define DEBUG_METIS false

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ███████ ████████ ██    ██ ██████  
//  ██      ██         ██    ██    ██ ██   ██ 
//  ███████ █████      ██    ██    ██ ██████  
//       ██ ██         ██    ██    ██ ██      
//  ███████ ███████    ██     ██████  ██      
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::setup( const NetworkParameters & network_parameters )
{
	//
	algorithm                  = network_parameters.algorithm ;
	objective                  = network_parameters.objective ;
	num_iter_save              = network_parameters.num_iter_save ;
	num_iter_print             = network_parameters.num_iter_print ;
	num_iter_max               = network_parameters.num_iter_max ;
	include_force_fixed_nodes  = network_parameters.include_force_fixed_nodes ;
	use_numerical_hessian      = network_parameters.use_numerical_hessian ;
	tolerance_change_objective = network_parameters.tolerance_change_objective ;
	tolerance_sum_net_force    = network_parameters.tolerance_sum_net_force ;
	//
	dir_input           = network_parameters.dir_input ;
	dir_output          = network_parameters.dir_output ;
	dir_output_iter     = dir_output + "iter" + FILESEP ;
	if( num_iter_save > 0 ) {
		make_dir( dir_output_iter ) ;
	}
	file_input_nodes    = dir_input  + "network_nodes.dat" ;
	file_input_springs  = dir_input  + "network_springs.dat" ;
	file_output_nodes   = dir_output + "network_nodes.dat" ;
	file_output_springs = dir_output + "network_springs.dat" ;
	std::random_device rd ;
	rng.seed( rd()  ) ;
	uni_0_1 = std::uniform_real_distribution<T>( 0.0,+1.0) ;
	uni_1_1 = std::uniform_real_distribution<T>(-1.0,+1.0) ;
	num_points  = network_parameters.num_points ;
	num_springs = network_parameters.num_springs ;
	points.resize( num_points ) ;
	springs.resize( num_springs ) ;
	nodes.resize( num_points ) ;
	load_network_binary( file_input_nodes.c_str() , file_input_springs.c_str() ) ;
	construct_network() ;
	setup_parallel( network_parameters ) ;
	partition_domain() ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██ 
// ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██            ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██  
// ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██            ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████   
// ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██            ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██  
//  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██    ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::construct_network( void )
{
	for( std::size_t p = 0 ; p < num_points ; ++p ) {
		nodes[p].point = &points[p] ;
		points[p].not_referenced = true ;
	}
	springs_used.clear() ;
	springs_used.reserve( num_springs ) ;
	for( std::size_t s = 0 ; s < num_springs ; ++s ) {
		std::size_t p_start = springs[s].start - &points[0] ;
		std::size_t p_end   = springs[s].end   - &points[0] ;
		// ignore link if spring has no force-length relationships, or if both points are fixed in all dimensions
		bool both_points_fixed = points[p_start].fixed_all_dim
		                      && points[p_end  ].fixed_all_dim ;
		bool none_force_length = springs[s].force_length_type_tension    ==Spring<T,N>::ForceLengthRelationship::none
		                      && springs[s].force_length_type_compression==Spring<T,N>::ForceLengthRelationship::none ;
		if( !both_points_fixed && !none_force_length ) {
			nodes[p_start].links.push_back( (Link){ -1 , &points[p_end  ] , &springs[s] , static_cast<T>(+1) } ) ;
			nodes[p_end  ].links.push_back( (Link){ -1 , &points[p_start] , &springs[s] , static_cast<T>(-1) } ) ;
			points[p_start].not_referenced = false ;
			points[p_end  ].not_referenced = false ;
			springs_used.push_back( &springs[s] ) ;
		}
	}
	nodes.erase( std::remove_if( nodes.begin() , nodes.end() , []( Node n ){ return n.point->fixed_all_dim  ; } ) , nodes.end() ) ;
	nodes.erase( std::remove_if( nodes.begin() , nodes.end() , []( Node n ){ return n.point->not_referenced ; } ) , nodes.end() ) ;
	for( std::size_t i = 0 ; i < nodes.size() ; ++i ) {
		nodes[i].node_index = i ;
		nodes[i].point->node_index = i ;
	}
	iterNode n ;
	iterNode n_end ;
	iterLink l ;
	iterLink l_end ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			l->node_index = l->point->node_index ;
		}
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ███████ ████████ ██    ██ ██████          ██████   █████  ██████   █████  ██      ██      ███████ ██      
//  ██      ██         ██    ██    ██ ██   ██         ██   ██ ██   ██ ██   ██ ██   ██ ██      ██      ██      ██      
//  ███████ █████      ██    ██    ██ ██████          ██████  ███████ ██████  ███████ ██      ██      █████   ██      
//       ██ ██         ██    ██    ██ ██              ██      ██   ██ ██   ██ ██   ██ ██      ██      ██      ██      
//  ███████ ███████    ██     ██████  ██      ███████ ██      ██   ██ ██   ██ ██   ██ ███████ ███████ ███████ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::setup_parallel( const NetworkParameters & network_parameters )
{
	//
	if( network_parameters.parallelism_enabled ) {
		if( network_parameters.user_num_threads > 1 ) {
			if( num_springs > 0 ) { // CHANGE BACK TO 1e3 *** TEST
				num_threads = static_cast<int>(network_parameters.user_num_threads) ;
				std::cout << "Multi-threaded: " << num_threads << " threads enabled." << std::endl ;
			} else {
				num_threads = 1 ;
				std::cout << "Single-threaded: Network contains <= 10^4 springs." << std::endl ;
			}
		} else {
			num_threads = 1 ;
			std::cout << "Single-threaded execution requested." << std::endl ;
		}
	} else {
		num_threads = 1 ;
		std::cout << "Single-threaded: Parallelism not enabled." << std::endl ;
	}
	parallelism_enabled = num_threads > 1 ;
	omp_set_num_threads( num_threads ) ;
	if( parallelism_enabled ) {
		omp_set_dynamic(false) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████   █████  ██████  ████████ ██ ████████ ██  ██████  ███    ██         ██████   ██████  ███    ███  █████  ██ ███    ██ 
//  ██   ██ ██   ██ ██   ██    ██    ██    ██    ██ ██    ██ ████   ██         ██   ██ ██    ██ ████  ████ ██   ██ ██ ████   ██ 
//  ██████  ███████ ██████     ██    ██    ██    ██ ██    ██ ██ ██  ██         ██   ██ ██    ██ ██ ████ ██ ███████ ██ ██ ██  ██ 
//  ██      ██   ██ ██   ██    ██    ██    ██    ██ ██    ██ ██  ██ ██         ██   ██ ██    ██ ██  ██  ██ ██   ██ ██ ██  ██ ██ 
//  ██      ██   ██ ██   ██    ██    ██    ██    ██  ██████  ██   ████ ███████ ██████   ██████  ██      ██ ██   ██ ██ ██   ████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::partition_domain( void )
{
	// prepare for compressed sparse row (csr) format of node-spring graph adjacency
	// have to eliminate coincident springs (multiple springs between same nodes) before passing to METIS
	// to account for computational cost, add coincident springs to edge weight
	// add weight to node if it has springs connected to fixed nodes (-1 node index)
	std::vector< std::vector< std::pair<idx_t,idx_t> > > adjncy_tmp( nodes.size() ) ;
	std::vector< idx_t > vwgt_tmp( nodes.size() , 0 ) ;
	std::size_t num_adj = 0 ;
	for( std::size_t i = 0 ; i < nodes.size() ; ++i ) {
		iterNode n = nodes.begin() + i ;
		iterLink l ;
		iterLink l_end ;
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			if( l->node_index < 0 ) {
				++ vwgt_tmp[i] ;
			} else {
				std::size_t j ;
				for( j = 0 ; j < adjncy_tmp[i].size() ; ++j ){
					if( adjncy_tmp[i][j].first == l->node_index ) {
						++ adjncy_tmp[i][j].second ;
						break ;
					}
				}
				if( j == adjncy_tmp[i].size() ) {
					adjncy_tmp[i].push_back( std::make_pair( l->node_index , 1 ) ) ;
					++ num_adj ;
				}
			}
		}
		std::sort( adjncy_tmp[i].begin() , adjncy_tmp[i].end() ) ;
	}
	// domain/graph partitioning via METIS library
	idx_t nvtxs = adjncy_tmp.size() ;
	idx_t ncon = 1 ;
	std::vector<idx_t> xadj( nvtxs + 1 ) ;
	std::vector<idx_t> adjncy( num_adj ) ;
	std::vector<idx_t> vwgt( nvtxs * ncon ) ; // this does not seem to help.
	std::vector<idx_t> adjwgt( num_adj ) ;
	std::size_t k = 0 ;
	xadj[0] = 0 ;
	for( std::size_t i = 0 ; i < adjncy_tmp.size() ; ++i ) {
		xadj[i+1] = xadj[i] + adjncy_tmp[i].size() ;
		vwgt[i] = vwgt_tmp[i] ;
		for( std::size_t j = 0 ; j < adjncy_tmp[i].size() ; ++j ) {
			adjncy[k] = adjncy_tmp[i][j].first ;
			adjwgt[k] = adjncy_tmp[i][j].second ;
			k++ ;
		}
	}
	idx_t nparts = num_threads ;
	idx_t options [METIS_NOPTIONS] ;
	idx_t objval ;
	std::vector<idx_t> part( nvtxs ) ;
	METIS_SetDefaultOptions(options) ;
	options[METIS_OPTION_PTYPE    ] = METIS_PTYPE_KWAY ;
	options[METIS_OPTION_OBJTYPE  ] = METIS_OBJTYPE_CUT ; // only for Kway
	options[METIS_OPTION_CTYPE    ] = METIS_CTYPE_SHEM ;
	options[METIS_OPTION_IPTYPE   ] = METIS_IPTYPE_EDGE ;
	options[METIS_OPTION_RTYPE    ] = METIS_RTYPE_GREEDY ;
	options[METIS_OPTION_NCUTS    ] = 1 ; // only for Kway
	options[METIS_OPTION_NUMBERING] = 0 ;
	options[METIS_OPTION_NITER    ] = 100 ;
	options[METIS_OPTION_CONTIG   ] = 1 ; // only for Kway
	options[METIS_OPTION_COMPRESS ] = 0 ;
	options[METIS_OPTION_UFACTOR  ] = 1 ;
	if( DEBUG_METIS ) {
		options[METIS_OPTION_DBGLVL] = METIS_DBG_COARSEN | METIS_DBG_REFINE ;
	} else {
		options[METIS_OPTION_DBGLVL] = 0 ;
	}
	if( nparts > 1 ) {
		METIS_PartGraphKway( &nvtxs , &ncon , &xadj.front() , &adjncy.front() , NULL  , NULL , &adjwgt.front() , &nparts , NULL , NULL , options , &objval , &part.front() ) ;
	} else {
		for( idx_t i = 0 ; i < nvtxs ; ++i ) {
			part[i] = 0 ;
		}
	}
	if( DEBUG_METIS ) {
		std::cout << "\nxadj\n" ;
		for( idx_t i = 0 ; i < nvtxs ; ++i ) {
			std::cout << "  " << xadj[i] ;
			if( ((i+1)%100) == 0 ) {
				std::cout << std::endl ;
			}
		}
		std::cout << std::endl ;
		std::cout << "\nadjncy\n" ;
		for( std::size_t i = 0 ; i < num_adj ; ++i ) {
			std::cout << "  " << adjncy[i] ;
			if( ((i+1)%100) == 0 ) {
				std::cout << std::endl ;
			}
		}
		std::cout << std::endl ;
		std::cout << "\npart\n" ;
		for( idx_t i = 0 ; i < nvtxs ; ++i ) {
			std::cout << part[i] << "  " ;
			if( ((i+1)%100) == 0 ) {
				std::cout << std::endl ;
			}
		}
		std::cout << std::endl ;
	}
	// create referencing structure for parallel execution of partitions
	thread_points.resize( num_threads ) ;
	thread_springs.resize( num_threads ) ;
	for( idx_t i = 0 ; i < nvtxs ; ++i ) {
		iterNode n = nodes.begin() + i ;
		iterLink l ;
		iterLink l_end ;
		// change shared to true if any other partition needs its position data
		bool shared = false ;
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			std::size_t spring_index = l->spring - &springs[0] ;
			if( l->node_index < 0 ) {
				// other point is a fixed point
				// current point is the only node with this spring
				// current point's partition must be responsible for updates
				thread_springs[ part[i] ].push_back( std::make_pair( spring_index , true ) ) ;
			} else if( part[i] != part[l->node_index] ) {
				// other point is in a different partition
				// both partitions should include the spring
				// only one partition is responsible for updates
				shared = true ;
				bool responsible = ( l->spring_direction > 0 ) ;
				thread_springs[ part[i] ].push_back( std::make_pair( spring_index , responsible ) ) ;
			} else {
				// other point is in the same partition
				// only one point should include the spring
				// that point must be responsible for updates
				if( l->spring_direction > 0 ) {
					thread_springs[ part[i] ].push_back( std::make_pair( spring_index , true ) ) ;
				}
			}
		}
		thread_points[ part[i] ].push_back( std::make_pair( nodes[i].point - &points[0] , shared ) ) ;
	}
	// print important infor for load balancing
	std::cout << "num of springs / points per thread:\n" ;
	for( idx_t i = 0 ; i < nparts ; ++i ) {
		std::cout << thread_springs[i].size() << "  " << thread_points[i].size() << "\n" ;
	}
	std::cout << std::endl ;

	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ███████ ████████ ██    ██ ██████          ███████ ██   ██  █████  ██████  ███████ ██████          ██████   █████  ████████  █████  
//  ██      ██         ██    ██    ██ ██   ██         ██      ██   ██ ██   ██ ██   ██ ██      ██   ██         ██   ██ ██   ██    ██    ██   ██ 
//  ███████ █████      ██    ██    ██ ██████          ███████ ███████ ███████ ██████  █████   ██   ██         ██   ██ ███████    ██    ███████ 
//       ██ ██         ██    ██    ██ ██                   ██ ██   ██ ██   ██ ██   ██ ██      ██   ██         ██   ██ ██   ██    ██    ██   ██ 
//  ███████ ███████    ██     ██████  ██      ███████ ███████ ██   ██ ██   ██ ██   ██ ███████ ██████  ███████ ██████  ██   ██    ██    ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::setup_shared_data( void )
{
	// initialize with original shared location of all points
	points_shared_location.resize(num_points) ;
	for( std::size_t p = 0 ; p < points.size() ; ++p ) {
		points_shared_location[p] = &points[p] ;
	}
	// assign subsets of modifiable data shared across threads
	// change locations of thread-private points to shared location, or NULL if not shared
	points_shared.resize(num_threads) ;
	for( std::size_t thread_id = 0 ; thread_id < num_threads ; ++thread_id ) {
		points_shared[thread_id].reserve( thread_points[thread_id].size() ) ;
		for( std::size_t i = 0 ; i < thread_points[thread_id].size() ; ++i ) {
			std::size_t p = thread_points[thread_id][i].first ;
			if( thread_points[thread_id][i].second ) { // shared?
				points_shared[thread_id].push_back( std::make_pair( points[p] , i ) ) ;
			}
		}
	}
	for( std::size_t thread_id = 0 ; thread_id < num_threads ; ++thread_id ) {
		std::size_t j = 0 ;
		for( std::size_t i = 0 ; i < thread_points[thread_id].size() ; ++i ) {
			std::size_t p = thread_points[thread_id][i].first ;
			if( thread_points[thread_id][i].second ) { // shared?
				points_shared_location[p] = & points_shared[thread_id][j].first ; // avoid pointer invalidation!
				++j ;
			} else {
				points_shared_location[p] = NULL ;
			}
		}
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ███████ ████████ ██    ██ ██████          ██       ██████   ██████  █████  ██              ██████   █████  ████████  █████  
//  ██      ██         ██    ██    ██ ██   ██         ██      ██    ██ ██      ██   ██ ██              ██   ██ ██   ██    ██    ██   ██ 
//  ███████ █████      ██    ██    ██ ██████          ██      ██    ██ ██      ███████ ██              ██   ██ ███████    ██    ███████ 
//       ██ ██         ██    ██    ██ ██              ██      ██    ██ ██      ██   ██ ██              ██   ██ ██   ██    ██    ██   ██ 
//  ███████ ███████    ██     ██████  ██      ███████ ███████  ██████   ██████ ██   ██ ███████ ███████ ██████  ██   ██    ██    ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::setup_local_data( const int thread_id ,
										   std::vector<Point<T,N>> & points_local ,
										   std::vector<Spring<T,N>> & springs_local ,
										   std::vector< std::vector< Link > > & links_local ,
										   std::vector< std::pair< Point<T,N> * , Point<T,N> * > > & points_update )
{
	// track where local points get stored by original index, NULL if not local
	std::vector< Point<T,N> * > points_private_location( num_points , NULL ) ;
	// subset of modifiable points local/private to this thread
	points_local.resize( thread_points[thread_id].size() ) ; // avoid pointer invalidation! ensure capacity >= max size
	for( std::size_t i = 0 ; i < thread_points[thread_id].size() ; ++i ) {
		std::size_t p = thread_points[thread_id][i].first ;
		points_local[i] = points[p] ;
		points_private_location[p] = & points_local[i] ; // avoid pointer invalidation!
	}
	// local update list for points_local to points_shared
	points_update.resize( points_shared[thread_id].size() ) ;
	for( std::size_t i = 0 ; i < points_shared[thread_id].size() ; ++i ) {
		points_update[i] = std::make_pair( & points_local[ points_shared[thread_id][i].second ] , & points_shared[thread_id][i].first ) ;
	}
	// subset of modifiable springs local/private to this thread
	// reset local spring pointers for local points, fixed points, and shared points
	// also create a vector of local link indexes for each local point, for local hessian calculation
	springs_local.resize( thread_springs[thread_id].size() ) ;
	links_local.resize( points_local.size() ) ;
	for( std::size_t i = 0 ; i < thread_springs[thread_id].size() ; ++i ) {
		std::size_t s = thread_springs[thread_id][i].first ;
		springs_local[i] = springs[s] ;
		std::size_t i_start = springs[s].start - &points[0] ;
		std::size_t i_end   = springs[s].end   - &points[0] ;
		springs_local[i].start = ( points_private_location[ i_start ] != NULL ) ? points_private_location[ i_start ] : points_shared_location[ i_start ] ; // avoid pointer invalidation!
		springs_local[i].end   = ( points_private_location[ i_end   ] != NULL ) ? points_private_location[ i_end   ] : points_shared_location[ i_end   ] ; // avoid pointer invalidation!
		int p_start = ( points_private_location[ i_start ] != NULL ) ? points_private_location[i_start]-&points_local[0] : -1 ;
		int p_end   = ( points_private_location[ i_end   ] != NULL ) ? points_private_location[i_end  ]-&points_local[0] : -1 ;
		if( points_private_location[ i_start ] != NULL ) {
			if( points_private_location[ i_end   ] != NULL ) {
				links_local[p_start].push_back( (Link){ p_end   , points_private_location[i_end  ] , &springs_local[i] , static_cast<T>(+1) } ) ; // avoid pointer invalidation!
			} else {
				links_local[p_start].push_back( (Link){ p_end   ,  points_shared_location[i_end  ] , &springs_local[i] , static_cast<T>(+1) } ) ; // avoid pointer invalidation!
			}
		}
		if( points_private_location[ i_end   ] != NULL ) {
			if( points_private_location[ i_start ] != NULL ) {
				links_local[p_end  ].push_back( (Link){ p_start , points_private_location[i_start] , &springs_local[i] , static_cast<T>(-1) } ) ; // avoid pointer invalidation!
			} else {
				links_local[p_end  ].push_back( (Link){ p_start ,  points_shared_location[i_start] , &springs_local[i] , static_cast<T>(-1) } ) ; // avoid pointer invalidation!
			}
		}
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██    ██ ██████  ██████   █████  ████████ ███████         ██████   ██████  ██ ███    ██ ████████ ███████         ███████ ██   ██  █████  ██████  ███████ ██████  
//  ██    ██ ██   ██ ██   ██ ██   ██    ██    ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ██      ██   ██ ██   ██ ██   ██ ██      ██   ██ 
//  ██    ██ ██████  ██   ██ ███████    ██    █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         ███████ ███████ ███████ ██████  █████   ██   ██ 
//  ██    ██ ██      ██   ██ ██   ██    ██    ██              ██      ██    ██ ██ ██  ██ ██    ██         ██              ██ ██   ██ ██   ██ ██   ██ ██      ██   ██ 
//   ██████  ██      ██████  ██   ██    ██    ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ███████ ██   ██ ██   ██ ██   ██ ███████ ██████  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::update_points_shared( const std::vector< std::pair< Point<T,N> * , Point<T,N> * > > & points_update )
{
	// synchronize shared point positions
	#pragma omp critical
	for( std::size_t i = 0 ; i < points_update.size() ; ++i ) {
		points_update[i].second->position = points_update[i].first->position ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██    ██ ██████  ██████   █████  ████████ ███████         ██████   ██████  ██ ███    ██ ████████ ███████          █████  ██      ██      
//  ██    ██ ██   ██ ██   ██ ██   ██    ██    ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ██   ██ ██      ██      
//  ██    ██ ██████  ██   ██ ███████    ██    █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         ███████ ██      ██      
//  ██    ██ ██      ██   ██ ██   ██    ██    ██              ██      ██    ██ ██ ██  ██ ██    ██         ██         ██   ██ ██      ██      
//   ██████  ██      ██████  ██   ██    ██    ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ██   ██ ███████ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::update_points_all( const std::vector< Point<T,N> > & points_subset )
{
	// synchronize all point positions
	for( std::size_t i = 0 ; i < points_subset.size() ; ++i ) {
		nodes[ points_subset[i].node_index ].point->position = points_subset[i].position ;
	}
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::update_points_all( const std::vector< Point<T,N> * > & points_subset )
{
	// synchronize all point positions
	for( std::size_t i = 0 ; i < points_subset.size() ; ++i ) {
		nodes[ points_subset[i]->node_index ].point->position = points_subset[i]->position ;
	}
	return ;
}
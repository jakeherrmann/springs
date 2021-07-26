#include "springs.hpp"
#include "vectors_nd.hpp"
#include "spmat.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>

#include <omp.h>
#include <stdlib.h>

#ifdef _WIN32
	#define FILESEP '\\'
#else
	#define FILESEP '/'
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███    ███ ██ ███    ██ ██ ███    ███ ██ ███████ ███████         ███████ ███    ██ ███████ ██████   ██████  ██    ██         ███    ██ ███████ ██     ██ ████████  ██████  ███    ██ 
//  ████  ████ ██ ████   ██ ██ ████  ████ ██    ███  ██              ██      ████   ██ ██      ██   ██ ██        ██  ██          ████   ██ ██      ██     ██    ██    ██    ██ ████   ██ 
//  ██ ████ ██ ██ ██ ██  ██ ██ ██ ████ ██ ██   ███   █████           █████   ██ ██  ██ █████   ██████  ██   ███   ████           ██ ██  ██ █████   ██  █  ██    ██    ██    ██ ██ ██  ██ 
//  ██  ██  ██ ██ ██  ██ ██ ██ ██  ██  ██ ██  ███    ██              ██      ██  ██ ██ ██      ██   ██ ██    ██    ██            ██  ██ ██ ██      ██ ███ ██    ██    ██    ██ ██  ██ ██ 
//  ██      ██ ██ ██   ████ ██ ██      ██ ██ ███████ ███████ ███████ ███████ ██   ████ ███████ ██   ██  ██████     ██    ███████ ██   ████ ███████  ███ ███     ██     ██████  ██   ████ 
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::minimize_energy_newton( void )
{
	// copy current state to initial, previous, and best states
	update_springs( springs_used ) ;
	update_forces() ;
	energy = total_energy() ;
	get_net_force_mag() ;
	T obj = get_objective() ;
	T obj_init = obj ;
	T obj_prev = obj ;
	//
	std::size_t num_iter_zero_change = 0 ;
	std::size_t num_iter_zero_change_max = 3 ;
	T step_size_reduction = 0.9 ; // in range (0,1)
	T step_size_min = 1E-16 ;
	T step_size_max = 0.10 ;
	T change_obj ;
	T mean_obj_change = static_cast<T>(0.0) ;
	T step_size_shared_max ;
	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 500 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :   0 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :   1 ;
	//
	std::size_t iter_parallel_max = 100 ;
	//
	std::size_t iter = 0 ;
	//
	T local_tolerance_sum_net_force    = ( tolerance_sum_net_force    > 0.0 ) ? tolerance_sum_net_force    : 1E-12 ;
	T local_tolerance_change_objective = ( tolerance_change_objective > 0.0 ) ? tolerance_change_objective : 1E-12 ;
	bool break_num_iter = false ;
	bool break_change_objective = false ;
	bool break_step_size = false ;
	bool break_sum_net_force_magnitude = false ;
	bool break_zero_change = false ;
	bool break_isnan_obj = false ;
	bool break_any = false ;

	//
	std::cout
		<< "  " << std::setw(10) << "E_curr"
		<< "  " << std::setw(10) << "E_prev"
		<< "  " << std::setw(10) << "E_change"
		<< "  " << std::setw(10) << "stepSize"
		<< "  " << std::setw(10) << "sum_F_net"
		<< "  " << std::setw(10) << "time"
		<< std::endl ;
	//
	double time_start = omp_get_wtime() ;
	double time_prev = time_start ;
	double time_curr ;

	// assign subsets of modifiable points shared by partitions
	setup_shared_data() ;
	T step_size_shared              [ num_threads ] ;
	T total_energy_shared           [ num_threads ] ;
	T max_net_force_magnitude_shared[ num_threads ] ;
	T sum_net_force_magnitude_shared[ num_threads ] ;
	#pragma omp parallel
	{
		// assign subsets of modifiable data private to each thread
		const int thread_id = omp_get_thread_num() ;
		std::vector< Point<T,N> > points_local ;
		std::vector< Spring<T,N> > springs_local ;
		std::vector< std::vector< Link > > links_local ;
		std::vector< std::pair< Point<T,N> * , Point<T,N> * > > points_update ;
		setup_local_data( thread_id , points_local , springs_local , links_local , points_update ) ;
		std::vector< Point<T,N> > points_local_prev = points_local ;
		std::default_random_engine rng_local{ std::random_device{}() } ;
		KleinSummer<T> ksum_local ;
		std::vector<T> step_direction_local ;
		spmat<T> hessian ;
		T total_energy_local ;
		T max_net_force_magnitude_local ;
		T sum_net_force_magnitude_local ;
		T step_size_local ;
		// begin newton method descent
		while( !break_any ) {
			#pragma omp barrier
			get_net_force_mag( points_local , ksum_local , max_net_force_magnitude_local , sum_net_force_magnitude_local ) ;
			total_energy_local = total_energy( springs_local , thread_springs[thread_id] , ksum_local ) ;
			T obj_local = get_objective( total_energy_local , sum_net_force_magnitude_local , max_net_force_magnitude_local ) ;
			for( std::size_t iter_parallel = 0 ; iter_parallel < iter_parallel_max ; ++iter_parallel ) {
				// basic line search
				T obj_prev_local = obj_local ;
				T step_size_local = 0.05 / (step_size_reduction*step_size_reduction) ;
				points_local_prev = points_local ;
				compute_newton_step_direction( points_local , springs_local , links_local , step_direction_local ) ;
				do {
					step_size_local *= step_size_reduction ;
					step_size = ( step_size > step_size_max ) ? step_size_max : step_size ;
					points_local = points_local_prev ;
					move_points_newton( step_size_local , points_local , step_direction_local ) ;
					update_springs( springs_local ) ;
					update_forces( points_local , springs_local ) ;
					get_net_force_mag( points_local , ksum_local , max_net_force_magnitude_local , sum_net_force_magnitude_local ) ;
					total_energy_local = total_energy( springs_local , thread_springs[thread_id] , ksum_local ) ;
					obj_local = get_objective( total_energy_local , sum_net_force_magnitude_local , max_net_force_magnitude_local ) ;
				} while( (obj_local>=obj_prev_local) && (step_size_local>step_size_min) ) ;
				update_points_shared( points_update ) ;
				update_springs( springs_local ) ;
				update_forces( points_local , springs_local ) ;
				get_net_force_mag( points_local , ksum_local , max_net_force_magnitude_local , sum_net_force_magnitude_local ) ;
				total_energy_local = total_energy( springs_local , thread_springs[thread_id] , ksum_local ) ;
			}
			step_size_shared              [ thread_id ] = step_size_local               ;
			total_energy_shared           [ thread_id ] = total_energy_local            ;
			max_net_force_magnitude_shared[ thread_id ] = max_net_force_magnitude_local ;
			sum_net_force_magnitude_shared[ thread_id ] = sum_net_force_magnitude_local ;
			#pragma omp barrier
			#pragma omp single
			{
				// go back and allocate partition-specific shared objectives below
				// merge partition-specific objectives into one total object
				energy = std::accumulate( &total_energy_shared[0] ,
										  &total_energy_shared[0] + num_threads ,
										  static_cast<T>(0) ) ;
				max_net_force_magnitude = *std::max_element( &max_net_force_magnitude_shared[0] ,
															 &max_net_force_magnitude_shared[0] + num_threads ) ;
				sum_net_force_magnitude = std::accumulate( &sum_net_force_magnitude_shared[0] ,
														   &sum_net_force_magnitude_shared[0] + num_threads ,
														   static_cast<T>(0) ) ;
				step_size_shared_max = *std::max_element( &step_size_shared[0] ,
														  &step_size_shared[0] + num_threads ) ;
				// evaluate objective function at new configuration
				obj_prev = obj ;
				obj = get_objective() ;
				mean_obj_change += (obj-obj_prev) ;
				//
				if( (local_num_iter_print>0) && (iter%local_num_iter_print==0) ) {
					mean_obj_change /= static_cast<T>(local_num_iter_print) ;
					time_curr = omp_get_wtime() ;
					std::cout
						<< std::setprecision(3)
						<< std::scientific
						<< "  " << std::setw(10) << obj
						<< "  " << std::setw(10) << obj_prev
						<< "  " << std::setw(10) << mean_obj_change
						<< "  " << std::setw(10) << step_size_shared_max
						<< "  " << std::setw(10) << sum_net_force_magnitude
						<< "  " << std::setw(10) << time_curr - time_prev
						<< std::endl ;
					time_prev = time_curr ;
					mean_obj_change = static_cast<T>(0.0) ;
				}
				// MUST SYNCHRONIZE ALL POINT POSITIONS BEFORE SAVING!
				/*
				if( (local_num_iter_save>0) && (iter%local_num_iter_save==0) ) {
					std::string dir_output_iter_curr = dir_output_iter + "iter_" + std::to_string(iter) + FILESEP ;
					make_dir( dir_output_iter_curr ) ;
					save_network_binary(
						(dir_output_iter_curr+"network_nodes.dat").c_str()   ,
						(dir_output_iter_curr+"network_springs.dat").c_str() ) ;
				}
				//*/
				//
				change_obj = (obj_prev-obj) / (obj_init-obj_prev) ;
				if( obj == obj_prev ) {
					++num_iter_zero_change ;
				}
				// stopping conditions
				++iter ;
				break_num_iter = iter >= local_num_iter_max ;
				break_change_objective = change_obj < local_tolerance_change_objective ;
				break_step_size = step_size_shared_max <= step_size_min ;
				break_sum_net_force_magnitude = sum_net_force_magnitude < local_tolerance_sum_net_force ;
				break_zero_change = num_iter_zero_change >= num_iter_zero_change_max ;
				break_isnan_obj = std::isnan( obj ) ;
				break_any = break_num_iter ||
					break_change_objective ||
					break_step_size ||
					break_sum_net_force_magnitude ||
					break_zero_change ||
					break_isnan_obj ;
			}
			//
			if( break_any ) {
				// synchronize best point positions to shared data
				#pragma omp critical
				{
					update_points_all( points_local ) ;
				}
				break ;
			}
		}
	}
	std::cout
		<< std::setprecision(3)
		<< std::scientific
		<< "  " << std::setw(10) << obj
		<< "  " << std::setw(10) << obj_prev
		<< "  " << std::setw(10) << obj-obj_prev
		<< "  " << std::setw(10) << step_size_shared_max
		<< "  " << std::setw(10) << sum_net_force_magnitude
		<< "  " << std::setw(10) << time_curr - time_prev
		<< std::endl ;

	// final update in best configuration
	update_springs() ;
	update_forces() ;
	energy = total_energy() ;
	get_net_force_mag() ;
	obj = get_objective() ;

	std::cout
		<< std::setprecision(3)
		<< std::scientific
		<< "  " << std::setw(10) << obj
		<< "  " << std::setw(10) << obj_prev
		<< "  " << std::setw(10) << obj-obj_prev
		<< "  " << std::setw(10) << step_size_shared_max
		<< "  " << std::setw(10) << sum_net_force_magnitude
		<< "  " << std::setw(10) << time_curr - time_prev
		<< std::endl ;

	std::cout << "solve time: " << omp_get_wtime() - time_start << std::endl ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███    ███  ██████  ██    ██ ███████         ██████   ██████  ██ ███    ██ ████████ ███████         ███    ██ ███████ ██     ██ ████████  ██████  ███    ██ 
// ████  ████ ██    ██ ██    ██ ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ████   ██ ██      ██     ██    ██    ██    ██ ████   ██ 
// ██ ████ ██ ██    ██ ██    ██ █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         ██ ██  ██ █████   ██  █  ██    ██    ██    ██ ██ ██  ██ 
// ██  ██  ██ ██    ██  ██  ██  ██              ██      ██    ██ ██ ██  ██ ██    ██         ██         ██  ██ ██ ██      ██ ███ ██    ██    ██    ██ ██  ██ ██ 
// ██      ██  ██████    ████   ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ██   ████ ███████  ███ ███     ██     ██████  ██   ████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_newton( const T & step_size )
{
	// move small displacement towards equilibrating position by 2nd order newton method
	/*
	Vector<T,N> small_displacement ;
	std::size_t d ;
	#pragma omp for schedule(static) private(small_displacement, d)
	for( std::size_t n = 0 ; n < nodes.size() ; ++n ) {
		for( d = 0 ; d < N ; ++d ) {
			small_displacement[d] = step_direction[ (nodes[n].node_index)*N + d ] ;
		}
		small_displacement *= step_size ;
		nodes[n].point->move( small_displacement ) ;
	}
	//*/
	//*
	Vector<T,N> small_displacement ;
	iterNode n ;
	iterNode n_end ;
	std::size_t d ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( d = 0 ; d < N ; ++d ) {
			small_displacement[d] = step_direction[ (n->node_index)*N + d ] ;
		}
		small_displacement *= step_size ;
		n->point->move( small_displacement ) ;
	}
	//*/
	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_newton( const T & step_size , std::vector< Point<T,N> > & points_subset , const std::vector<T> & step_direction )
{
	Vector<T,N> small_displacement ;
	std::size_t d ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		for( d = 0 ; d < N ; ++d ) {
			small_displacement[d] = step_direction[ p*N + d ] ;
		}
		small_displacement *= step_size ;
		points_subset[p].move( small_displacement ) ;
	}
	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_newton( const T & step_size , std::vector< Point<T,N> * > & points_subset , const std::vector<T> & step_direction )
{
	Vector<T,N> small_displacement ;
	std::size_t d ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		for( d = 0 ; d < N ; ++d ) {
			small_displacement[d] = step_direction[ p*N + d ] ;
		}
		small_displacement *= step_size ;
		points_subset[p]->move( small_displacement ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████          ██████  ██████   █████  ██████  ██ ███████ ███    ██ ████████ 
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ██       ██   ██ ██   ██ ██   ██ ██ ██      ████   ██    ██    
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ██   ███ ██████  ███████ ██   ██ ██ █████   ██ ██  ██    ██    
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██    ██ ██   ██ ██   ██ ██   ██ ██ ██      ██  ██ ██    ██    
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████  ██████  ██   ██ ██   ██ ██████  ██ ███████ ██   ████    ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_gradient( void )
{
	// evalute net force at each node (gradient of total spring network energy)
	std::size_t num_node = nodes.size() ;
	neg_gradient.resize( N*num_node ) ;
	iterNode n ;
	iterNode n_end ;
	std::size_t d ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( d = 0 ; d < N ; ++d ) {
			neg_gradient[ (n->node_index)*N + d ] = n->point->net_force[d] ;
		}
	}
	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_gradient( const std::vector< Point<T,N> > & points_subset , std::vector<T> & neg_gradient )
{
	// evalute net force at each node (gradient of total spring network energy)
	std::size_t num_node = points_subset.size() ;
	neg_gradient.resize( N*num_node ) ;
	std::size_t d ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		for( d = 0 ; d < N ; ++d ) {
			neg_gradient[ p*N + d ] = points_subset[p].net_force[d] ;
		}
	}
	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_gradient( const std::vector< Point<T,N> * > & points_subset , std::vector<T> & neg_gradient )
{
	// evalute net force at each node (gradient of total spring network energy)
	std::size_t num_node = points_subset.size() ;
	neg_gradient.resize( N*num_node ) ;
	std::size_t d ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		for( d = 0 ; d < N ; ++d ) {
			neg_gradient[ p*N + d ] = points_subset[p]->net_force[d] ;
		}
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████         ██   ██ ███████ ███████ ███████ ██  █████  ███    ██         ███    ██ ██    ██ ███    ███ ███████ ██████  ██  ██████  █████  ██      
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ██   ██ ██      ██      ██      ██ ██   ██ ████   ██         ████   ██ ██    ██ ████  ████ ██      ██   ██ ██ ██      ██   ██ ██      
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ███████ █████   ███████ ███████ ██ ███████ ██ ██  ██         ██ ██  ██ ██    ██ ██ ████ ██ █████   ██████  ██ ██      ███████ ██      
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██   ██ ██           ██      ██ ██ ██   ██ ██  ██ ██         ██  ██ ██ ██    ██ ██  ██  ██ ██      ██   ██ ██ ██      ██   ██ ██      
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████ ██   ██ ███████ ███████ ███████ ██ ██   ██ ██   ████ ███████ ██   ████  ██████  ██      ██ ███████ ██   ██ ██  ██████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_hessian_numerical( void )
{
	// evaluate spatial derivatives in nodal forces
	// using central difference, small step +/- in each direction
	// select only springs connected to node n to re-evaluate forces
	T delta_position = static_cast<T>(1.0e-4) ; // assuming network is rescaled
	T delta_position_half = delta_position * static_cast<T>(0.5) ;
	std::size_t num_node = nodes.size() ;
	std::vector<std::size_t> sp_row ;
	std::vector<std::size_t> sp_col ;
	std::vector<T> sp_val ;
	iterNode n ;
	iterNode n_end ;
	iterLink l ;
	iterLink l_end ;
	std::size_t d_p ;
	std::size_t d_f ;
	Vector<T,N> force_step_pos ;
	Vector<T,N> force_step_neg ;
	Vector<T,N> force_delta ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		// change in net force at n'th node w.r.t. small changes in n'th node position
		for( d_p = 0 ; d_p < N ; ++d_p ) {
			if( ! n->point->fixed_dim[d_p] ) {
				n->point->position[d_p] += delta_position_half ;
				force_step_pos = static_cast<T>(0.0) ;
				for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
					l->spring->spring_energy() ;
					if( l->spring_direction > static_cast<T>(0) ) {
						force_step_pos += l->spring->force ;
					} else {
						force_step_pos -= l->spring->force ;
					}
				}
				n->point->position[d_p] -= delta_position ;
				force_step_neg = static_cast<T>(0.0) ;
				for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
					l->spring->spring_energy() ;
					if( l->spring_direction > static_cast<T>(0) ) {
						force_step_neg += l->spring->force ;
					} else {
						force_step_neg -= l->spring->force ;
					}
				}
				n->point->position[d_p] += delta_position_half ;
				force_delta = force_step_pos - force_step_neg ;
				for( d_f = 0 ; d_f < N ; ++d_f ) {
					sp_row.push_back( (n->node_index)*N + d_f ) ;
					sp_col.push_back( (n->node_index)*N + d_p ) ;
					sp_val.push_back( force_delta[d_f] / delta_position ) ;
				}
			}
		}
		// change in net force at n'th node w.r.t. small changes in linked node positions
		// only consider moveable nodes (i.e., not fixed points)
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			if( l->node_index >= 0 ) {
				for( d_p = 0 ; d_p < N ; ++d_p ) {
					if( ! l->point->fixed_dim[d_p] ) {
						l->point->position[d_p] += delta_position_half ;
						l->spring->spring_energy() ;
						force_step_pos = l->spring->force ;
						l->point->position[d_p] -= delta_position ;
						l->spring->spring_energy() ;
						force_step_neg = l->spring->force ;
						l->point->position[d_p] += delta_position_half ;
						if( l->spring_direction > static_cast<T>(0) ) {
							force_delta = force_step_pos - force_step_neg ;
						} else {
							force_delta = force_step_neg - force_step_pos ;
						}
						for( d_f = 0 ; d_f < N ; ++d_f ) {
							sp_row.push_back( (n->node_index)*N + d_f ) ;
							sp_col.push_back( (l->node_index)*N + d_p ) ;
							sp_val.push_back( force_delta[d_f] / delta_position ) ;
						}
					}
				}	
			}
		}
	}
	hessian = typename spmat<T>::spmat( N*num_node , N*num_node ) ;
	hessian.set_val( sp_row , sp_col , sp_val ) ;
	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_hessian_numerical( std::vector< Point<T,N> > & points_subset , const std::vector<std::vector<Link>> & links_subset , spmat<T> & hessian_subset )
{
	// evaluate spatial derivatives in nodal forces
	// using central difference, small step +/- in each direction
	// select only springs connected to node n to re-evaluate forces
	T delta_position = static_cast<T>(1.0e-4) ; // assuming network is rescaled
	T delta_position_half = delta_position * static_cast<T>(0.5) ;
	std::size_t num_node = points_subset.size() ; // safe because fixed points are not assigned to subset partitions
	std::vector<std::size_t> sp_row ;
	std::vector<std::size_t> sp_col ;
	std::vector<T> sp_val ;
	iterLink l ;
	iterLink l_end ;
	std::size_t d_p ;
	std::size_t d_f ;
	Vector<T,N> force_step_pos ;
	Vector<T,N> force_step_neg ;
	Vector<T,N> force_delta ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		// change in net force at n'th node w.r.t. small changes in n'th node position
		for( d_p = 0 ; d_p < N ; ++d_p ) {
			if( ! points_subset[p].fixed_dim[d_p] ) {
				points_subset[p].position[d_p] += delta_position_half ;
				force_step_pos = static_cast<T>(0.0) ;
				for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
					l->spring->spring_energy() ;
					if( l->spring_direction > static_cast<T>(0) ) {
						force_step_pos += l->spring->force ;
					} else {
						force_step_pos -= l->spring->force ;
					}
				}
				points_subset[p].position[d_p] -= delta_position ;
				force_step_neg = static_cast<T>(0.0) ;
				for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
					l->spring->spring_energy() ;
					if( l->spring_direction > static_cast<T>(0) ) {
						force_step_neg += l->spring->force ;
					} else {
						force_step_neg -= l->spring->force ;
					}
				}
				points_subset[p].position[d_p] += delta_position_half ;
				force_delta = force_step_pos - force_step_neg ;
				for( d_f = 0 ; d_f < N ; ++d_f ) {
					sp_row.push_back( p*N + d_f ) ;
					sp_col.push_back( p*N + d_p ) ;
					sp_val.push_back( force_delta[d_f] / delta_position ) ;
				}
			}
		}
		// change in net force at n'th node w.r.t. small changes in linked node positions
		// only consider moveable nodes (i.e., not fixed points)
		// do not attempt to move points not in local partition (i.e., not null pointer)
		for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
			if( l->node_index >= 0 ) {
				for( d_p = 0 ; d_p < N ; ++d_p ) {
					if( l->point != NULL ) {
						if( ! l->point->fixed_dim[d_p] ) {
							l->point->position[d_p] += delta_position_half ;
							l->spring->spring_energy() ;
							force_step_pos = l->spring->force ;
							l->point->position[d_p] -= delta_position ;
							l->spring->spring_energy() ;
							force_step_neg = l->spring->force ;
							l->point->position[d_p] += delta_position_half ;
							if( l->spring_direction > static_cast<T>(0) ) {
								force_delta = force_step_pos - force_step_neg ;
							} else {
								force_delta = force_step_neg - force_step_pos ;
							}
							for( d_f = 0 ; d_f < N ; ++d_f ) {
								sp_row.push_back( p*N + d_f ) ;
								sp_col.push_back( (l->node_index)*N + d_p ) ;
								sp_val.push_back( force_delta[d_f] / delta_position ) ;
							}
						}
					}
				}	
			}
		}
	}
	hessian_subset = typename spmat<T>::spmat( N*num_node , N*num_node ) ;
	hessian_subset.set_val( sp_row , sp_col , sp_val ) ;
	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_hessian_numerical( const std::vector< Point<T,N> * > & points_subset , const std::vector<std::vector<Link>> & links_subset , spmat<T> & hessian_subset )
{
	// evaluate spatial derivatives in nodal forces
	// using central difference, small step +/- in each direction
	// select only springs connected to node n to re-evaluate forces
	T delta_position = static_cast<T>(1.0e-4) ; // assuming network is rescaled
	T delta_position_half = delta_position * static_cast<T>(0.5) ;
	std::size_t num_node = points_subset.size() ; // safe because fixed points are not assigned to subset partitions
	std::vector<std::size_t> sp_row ;
	std::vector<std::size_t> sp_col ;
	std::vector<T> sp_val ;
	iterLink l ;
	iterLink l_end ;
	std::size_t d_p ;
	std::size_t d_f ;
	Vector<T,N> force_step_pos ;
	Vector<T,N> force_step_neg ;
	Vector<T,N> force_delta ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		// change in net force at n'th node w.r.t. small changes in n'th node position
		for( d_p = 0 ; d_p < N ; ++d_p ) {
			if( ! points_subset[p]->fixed_dim[d_p] ) {
				points_subset[p]->position[d_p] += delta_position_half ;
				force_step_pos = static_cast<T>(0.0) ;
				for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) { ///
					l->spring->spring_energy() ;
					if( l->spring_direction > static_cast<T>(0) ) {
						force_step_pos += l->spring->force ;
					} else {
						force_step_pos -= l->spring->force ;
					}
				}
				points_subset[p]->position[d_p] -= delta_position ;
				force_step_neg = static_cast<T>(0.0) ;
				for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
					l->spring->spring_energy() ;
					if( l->spring_direction > static_cast<T>(0) ) {
						force_step_neg += l->spring->force ;
					} else {
						force_step_neg -= l->spring->force ;
					}
				}
				points_subset[p]->position[d_p] += delta_position_half ;
				force_delta = force_step_pos - force_step_neg ;
				for( d_f = 0 ; d_f < N ; ++d_f ) {
					sp_row.push_back( p*N + d_f ) ;
					sp_col.push_back( p*N + d_p ) ;
					sp_val.push_back( force_delta[d_f] / delta_position ) ;
				}
			}
		}
		// change in net force at n'th node w.r.t. small changes in linked node positions
		// only consider moveable nodes (i.e., not fixed points)
		// do not attempt to move points not in local partition (i.e., not null pointer)
		for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
			if( l->node_index >= 0 ) {
				for( d_p = 0 ; d_p < N ; ++d_p ) {
					if( l->point != NULL ) {
						if( ! l->point->fixed_dim[d_p] ) {
							l->point->position[d_p] += delta_position_half ;
							l->spring->spring_energy() ;
							force_step_pos = l->spring->force ;
							l->point->position[d_p] -= delta_position ;
							l->spring->spring_energy() ;
							force_step_neg = l->spring->force ;
							l->point->position[d_p] += delta_position_half ;
							if( l->spring_direction > static_cast<T>(0) ) {
								force_delta = force_step_pos - force_step_neg ;
							} else {
								force_delta = force_step_neg - force_step_pos ;
							}
							for( d_f = 0 ; d_f < N ; ++d_f ) {
								sp_row.push_back( p*N + d_f ) ;
								sp_col.push_back( (l->node_index)*N + d_p ) ;
								sp_val.push_back( force_delta[d_f] / delta_position ) ;
							}
						}
					}
				}	
			}
		}
	}
	hessian_subset = typename spmat<T>::spmat( N*num_node , N*num_node ) ;
	hessian_subset.set_val( sp_row , sp_col , sp_val ) ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████         ██   ██ ███████ ███████ ███████ ██  █████  ███    ██          █████  ███    ██  █████  ██      ██    ██ ████████ ██  ██████  █████  ██      
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ██   ██ ██      ██      ██      ██ ██   ██ ████   ██         ██   ██ ████   ██ ██   ██ ██       ██  ██     ██    ██ ██      ██   ██ ██      
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ███████ █████   ███████ ███████ ██ ███████ ██ ██  ██         ███████ ██ ██  ██ ███████ ██        ████      ██    ██ ██      ███████ ██      
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██   ██ ██           ██      ██ ██ ██   ██ ██  ██ ██         ██   ██ ██  ██ ██ ██   ██ ██         ██       ██    ██ ██      ██   ██ ██      
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████ ██   ██ ███████ ███████ ███████ ██ ██   ██ ██   ████ ███████ ██   ██ ██   ████ ██   ██ ███████    ██       ██    ██  ██████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_hessian_analytical( void )
{
	//
	iterNode n ;
	iterNode n_end ;
	iterLink l ;
	iterLink l_end ;
	std::size_t num_node = nodes.size() ;
	std::vector<std::size_t> sp_row ;
	std::vector<std::size_t> sp_col ;
	std::vector<T> sp_val ;
	//
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		// case: different nodes, any directions
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			if( l->node_index >= 0 ) {
				T value_0 = l->spring->rest_length / l->spring->length ;
				T value_1 = value_0 - static_cast<T>(1.0) ;
				T value_2 = value_0 / ( l->spring->length * l->spring->length ) ;
				for( std::size_t d_p = 0 ; d_p < N ; ++d_p ) {
					for( std::size_t d_f = 0 ; d_f < N ; ++d_f ) {
						T value_3 = ( n->point->position[d_p] - l->point->position[d_p] ) ;
						T value_4 = ( n->point->position[d_f] - l->point->position[d_f] ) ;
						T value_5 = l->spring->effective_spring_constant * ( value_1 - (value_2*value_3*value_4) ) ;
						sp_row.push_back( (n->node_index)*N + d_f ) ;
						sp_col.push_back( (l->node_index)*N + d_p ) ;
						sp_val.push_back( value_5 ) ;
					}
				}
			}
		}
		// case: same node, different directions
		for( std::size_t d_p = 0 ; d_p < N ; ++d_p ) {
			for( std::size_t d_f = 0 ; d_f < N ; ++d_f ) {
				T value_0 = static_cast<T>(0.0) ;
				for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
					T value_1 = l->spring->rest_length / l->spring->length ;
					T value_2 = static_cast<T>(1.0) - value_1 ;
					T value_3 = value_1 / ( l->spring->length * l->spring->length ) ;
					T value_4 = ( n->point->position[d_p] - l->point->position[d_p] ) ;
					T value_5 = ( n->point->position[d_f] - l->point->position[d_f] ) ;
					value_0 += l->spring->effective_spring_constant * ( value_2 + value_3*value_4*value_5 ) ;
				}
				sp_row.push_back( (n->node_index)*N + d_f ) ;
				sp_col.push_back( (n->node_index)*N + d_p ) ;
				sp_val.push_back( value_0 ) ;
			}
		}
	}
	hessian = typename spmat<T>::spmat( N*num_node , N*num_node ) ;
	hessian.set_val( sp_row , sp_col , sp_val ) ;
	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_hessian_analytical( std::vector< Point<T,N> > & points_subset , std::vector< Spring<T,N> > & springs_subset , const std::vector<std::vector<Link>> & links_subset , spmat<T> & hessian_subset )
{
	//
	iterLink l ;
	iterLink l_end ;
	std::size_t num_node = points_subset.size() ; // safe because fixed points are not assigned to subset partitions
	std::size_t num_spring = springs_subset.size() ;
	std::vector<std::size_t> sp_row ;
	std::vector<std::size_t> sp_col ;
	std::vector<T> sp_val ;
	//
	std::vector<std::pair<T,T>> hessian_coef ;
	hessian_coef.resize( num_spring ) ;
	for( std::size_t s = 0 ; s < springs_subset.size() ; ++s ) {
		hessian_coef[s] = springs_subset[s].hessian_coefs() ;
	}
	//
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
			std::size_t s = l->spring - &springs_subset[0] ;
			for( std::size_t d_p = 0 ; d_p < N ; ++d_p ) {
				T delta_p = points_subset[p].position[d_p] - l->point->position[d_p] ;
				sp_row.push_back( p*N + d_p ) ;
				sp_col.push_back( (l->node_index)*N + d_p ) ;
				sp_val.push_back( -( hessian_coefs[s].first + (delta_p*delta_p*hessian_coefs[s].second) ) ) ;
				for( std::size_t d_o = 0 ; d_o < (N-1) ; ++d_o ) {
					std::size_t d_f = ( d_p + d_o ) % N ;
					T delta_f = points_subset[p].position[d_f] - l->point->position[d_f] ;
					T val = -( delta_p * delta_f * hessian_coefs[s].second )  ;
					sp_row.push_back( p*N + d_p ) ;
					sp_col.push_back( (l->node_index)*N + d_f) ;
					sp_val.push_back( val ) ;
					sp_row.push_back( (l->node_index)*N + d_f) ;
					sp_col.push_back( p*N + d_p ) ;
					sp_val.push_back( val ) ;
				}
			}
		}
		for( std::size_t d_p = 0 ; d_p < N ; ++d_p ) {
			T val = static_cast<T>(0) ;
			for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
				std::size_t s = l->spring - &springs_subset[0] ;
				T delta_p = points_subset[p].position[d_p] - l->point->position[d_p] ;
				val += hessian_coefs[s].first - (delta_p*delta_p*hessian_coefs[s].second) ;
			}
			sp_row.push_back( p*N + d_p ) ;
			sp_col.push_back( p*N + d_p ) ;
			sp_val.push_back( val ) ;
			val = static_cast<T>(0) ;
			for( std::size_t d_o = 0 ; d_o < (N-1) ; ++d_o ) {
				std::size_t d_f = ( d_p + d_o ) % N ;
				for( l = links_subset[p].begin() , l_end = links_subset[p].end() ; l != l_end ; ++l ) {
					std::size_t s = l->spring - &springs_subset[0] ;
					T delta_f = points_subset[p].position[d_f] - l->point->position[d_f] ;
					val += (delta_p*delta_f*hessian_coefs[s].second) ;
				}
				sp_row.push_back( p*N + d_p ) ;
				sp_col.push_back( p*N + d_f ) ;
				sp_val.push_back( val ) ;
				sp_row.push_back( p*N + d_f ) ;
				sp_col.push_back( p*N + d_p ) ;
				sp_val.push_back( val ) ;
			}
		}
	}
	hessian_subset = typename spmat<T>::spmat( N*num_node , N*num_node ) ;
	hessian_subset.set_val( sp_row , sp_col , sp_val ) ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████         ███    ██ ███████ ██     ██ ████████  ██████  ███    ██         ███████ ████████ ███████ ██████          ██████  ██ ██████  ███████  ██████ ████████ ██  ██████  ███    ██ 
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ████   ██ ██      ██     ██    ██    ██    ██ ████   ██         ██         ██    ██      ██   ██         ██   ██ ██ ██   ██ ██      ██         ██    ██ ██    ██ ████   ██ 
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ██ ██  ██ █████   ██  █  ██    ██    ██    ██ ██ ██  ██         ███████    ██    █████   ██████          ██   ██ ██ ██████  █████   ██         ██    ██ ██    ██ ██ ██  ██ 
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██  ██ ██ ██      ██ ███ ██    ██    ██    ██ ██  ██ ██              ██    ██    ██      ██              ██   ██ ██ ██   ██ ██      ██         ██    ██ ██    ██ ██  ██ ██ 
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████ ██   ████ ███████  ███ ███     ██     ██████  ██   ████ ███████ ███████    ██    ███████ ██      ███████ ██████  ██ ██   ██ ███████  ██████    ██    ██  ██████  ██   ████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_newton_step_direction( void )
{
	// prepare gradient
	compute_gradient() ;
	step_direction = neg_gradient ;

	// prepare hessian
	if( use_numerical_hessian ) {
		compute_hessian_numerical() ;
	} else {
		compute_hessian_analytical() ;
	}

	// if any zero on the diagonal, remove row and column
	// set all row and col values 0, set diagonal value 1, set source vector 0
	T small_number = static_cast<T>(1000.0) * std::numeric_limits<T>::epsilon() ;
	for( std::size_t r = 0 ; r < hessian.get_numRow() ; ++r ) {
		if( std::abs(hessian.get_val(r,r)) <= small_number ) {
			for( std::size_t c = 0 ; c < hessian.get_numCol() ; ++c ) {
				hessian.set_val(r,c,static_cast<T>(0.0)) ;
				hessian.set_val(c,r,static_cast<T>(0.0)) ;
			}
			hessian.set_val(r,r,static_cast<T>(1.0)) ;
			neg_gradient[r] = static_cast<T>(0.0) ;
			step_direction[r] = static_cast<T>(0.0) ;
		}
	}

	// solve using gauss-seidel iterations with successive over-relaxation
	hessian.solve_SOR(
		neg_gradient ,
		step_direction ,
		static_cast<T>(1.0e-6) , // error tolerance
		static_cast<std::size_t>(1000) , // max number of iterations
		true ) ; // dynamically modify relaxation parameter?

	/*
	std::cout << "__________________________________" << std::endl ;
	for( iterNode n = nodes.begin() ; n != nodes.end() ; ++n ) {
		for( std::size_t d = 0 ; d < N ; ++d ) {
			std::cout << '\t' << step_direction[ (n->node_index)*N + d ] ;
		}
		std::cout << '\n' ;
	}
	std::cout << std::endl ;
	//*/

	//*
	// ensure the newton step is aligned with gradient
	T dot_product = 0.0 ;
	for( std::size_t r = 0 ; r < hessian.get_numRow() ; ++r ) {
		dot_product += neg_gradient[r] * step_direction[r] ;
	}
	if( std::isnan(dot_product) ) {
		//std::cout << "bad newton step: NAN" << std::endl ;
		step_direction = neg_gradient ;
	} else if( dot_product <= 0.0 ) {
		//std::cout << "bad newton step: opposite to gradient" << std::endl ;
		step_direction = neg_gradient ;
	} else {
		//std::cout << "good" << std::endl ;
	}
	//*/

	return ;
}
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_newton_step_direction( std::vector< Point<T,N> > & points_subset , std::vector< Spring<T,N> > & springs_subset , const std::vector<std::vector<Link>> & links_subset , std::vector<T> & step_direction )
{
	// prepare gradient
	std::vector<T> neg_gradient ;
	compute_gradient( points_subset , neg_gradient ) ;
	step_direction = neg_gradient ;

	// prepare hessian
	spmat<T> hessian_subset ;
	if( use_numerical_hessian ) {
		compute_hessian_numerical( points_subset , links_subset , hessian_subset ) ;
	} else {
		compute_hessian_analytical( points_subset , springs_subset , links_subset , hessian_subset ) ;
	}

	// if any zero on the diagonal, remove row and column
	// set all row and col values 0, set diagonal value 1, set source vector 0
	T small_number = static_cast<T>(1000.0) * std::numeric_limits<T>::epsilon() ;
	for( std::size_t r = 0 ; r < hessian_subset.get_numRow() ; ++r ) {
		if( std::abs(hessian_subset.get_val(r,r)) <= small_number ) {
			for( std::size_t c = 0 ; c < hessian_subset.get_numCol() ; ++c ) {
				hessian_subset.set_val(r,c,static_cast<T>(0.0)) ;
				hessian_subset.set_val(c,r,static_cast<T>(0.0)) ;
			}
			hessian_subset.set_val(r,r,static_cast<T>(1.0)) ;
			neg_gradient[r] = static_cast<T>(0.0) ;
			step_direction[r] = static_cast<T>(0.0) ;
		}
	}

	// solve using gauss-seidel iterations with successive over-relaxation
	hessian_subset.solve_SOR(
		neg_gradient ,
		step_direction ,
		static_cast<T>(1.0e-6) , // error tolerance
		static_cast<std::size_t>(1000) , // max number of iterations
		true ) ; // dynamically modify relaxation parameter?

	//*
	// ensure the newton step is aligned with gradient
	T dot_product = 0.0 ;
	for( std::size_t r = 0 ; r < hessian_subset.get_numRow() ; ++r ) {
		dot_product += neg_gradient[r] * step_direction[r] ;
	}
	if( std::isnan(dot_product) ) {
		//std::cout << "bad newton step: NAN" << std::endl ;
		step_direction = neg_gradient ;
	} else if( dot_product <= 0.0 ) {
		//std::cout << "bad newton step: opposite to gradient" << std::endl ;
		step_direction = neg_gradient ;
	} else {
		//std::cout << "good" << std::endl ;
	}
	//*/

	return ;
}
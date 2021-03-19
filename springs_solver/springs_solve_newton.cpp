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
//  ███    ███ ██ ███    ██ ██ ███    ███ ██ ███████ ███████         ███████ ███    ██ ███████ ██████   ██████  ██    ██         ███    ██ ███████ ██     ██ ████████  ██████  ███    ██ 
//  ████  ████ ██ ████   ██ ██ ████  ████ ██    ███  ██              ██      ████   ██ ██      ██   ██ ██        ██  ██          ████   ██ ██      ██     ██    ██    ██    ██ ████   ██ 
//  ██ ████ ██ ██ ██ ██  ██ ██ ██ ████ ██ ██   ███   █████           █████   ██ ██  ██ █████   ██████  ██   ███   ████           ██ ██  ██ █████   ██  █  ██    ██    ██    ██ ██ ██  ██ 
//  ██  ██  ██ ██ ██  ██ ██ ██ ██  ██  ██ ██  ███    ██              ██      ██  ██ ██ ██      ██   ██ ██    ██    ██            ██  ██ ██ ██      ██ ███ ██    ██    ██    ██ ██  ██ ██ 
//  ██      ██ ██ ██   ████ ██ ██      ██ ██ ███████ ███████ ███████ ███████ ██   ████ ███████ ██   ██  ██████     ██    ███████ ██   ████ ███████  ███ ███     ██     ██████  ██   ████ 
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::minimize_energy_newton( void )
{
	// copy current state to initial, previous, and best states
	points_init = points ;
	points_prev = points ;
	update_springs() ;
	update_forces() ; // also computes sum_net_force_magnitude
	T energy = get_objective() ;
	T energy_init = energy ;
	T energy_prev = energy ;

	//
	std::size_t num_iter_zero_change = 0 ;
	std::size_t num_iter_zero_change_max = 3 ;
	T step_size = 0.05 ;
	T step_size_reduction = 0.9 ; // in range (0,1)
	T step_size_min = 1E-16 ;
	T step_size_max = 0.10 ;
	T change_energy ;

	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 500 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :   0 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :   1 ;
	T local_tolerance_sum_net_force    = ( tolerance_sum_net_force    > 0.0 ) ? tolerance_sum_net_force    : 1E-12 ;
	T local_tolerance_change_objective = ( tolerance_change_objective > 0.0 ) ? tolerance_change_objective : 1E-12 ;
	bool break_change_objective = false ;
	bool break_step_size = false ;
	bool break_sum_net_force_magnitude = false ;
	bool break_zero_change = false ;
	bool break_isnan_energy = false ;

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

	//
	for( std::size_t iter = 0 ; iter < local_num_iter_max ; ++iter ) {

		// basic line search
		energy_prev = energy ;
		points_prev = points ;
		step_size /= step_size_reduction ;
		step_size /= step_size_reduction ;
		compute_newton_step_direction() ;

		do {
			step_size *= step_size_reduction ;
			step_size = ( step_size > step_size_max ) ? step_size_max : step_size ;
			points = points_prev ;
			move_points_newton( step_size ) ;
			update_springs() ;
			update_forces() ; // also computes sum_net_force_magnitude
			energy = get_objective() ;
		} while( (energy>=energy_prev) && (step_size>step_size_min) ) ;

		//
		change_energy = (energy_prev-energy) / (energy_init-energy_prev) ;
		if( energy == energy_prev ) {
			++num_iter_zero_change ;
		}

		if( (local_num_iter_print>0) && (iter%local_num_iter_print==0) ) {
			time_curr = omp_get_wtime() ;
			std::cout
				<< std::setprecision(3)
				<< std::scientific
				<< "  " << std::setw(10) << energy
				<< "  " << std::setw(10) << energy_prev
				<< "  " << std::setw(10) << energy-energy_prev
				<< "  " << std::setw(10) << step_size
				<< "  " << std::setw(10) << sum_net_force_magnitude
				<< "  " << std::setw(10) << time_curr - time_prev
				<< std::endl ;
			time_prev = time_curr ;
		}

		if( (local_num_iter_save>0) && (iter%local_num_iter_save==0) ) {
			std::string dir_output_iter_curr = dir_output_iter + "iter_" + std::to_string(iter) + FILESEP ;
			make_dir( dir_output_iter_curr ) ;
			save_network_binary(
				(dir_output_iter_curr+"network_nodes.dat").c_str()   ,
				(dir_output_iter_curr+"network_springs.dat").c_str() ) ;
		}

		// stopping conditions
		break_change_objective = change_energy < local_tolerance_change_objective ;\
		break_step_size = step_size <= step_size_min ;
		break_sum_net_force_magnitude = sum_net_force_magnitude < local_tolerance_sum_net_force ;
		break_zero_change = num_iter_zero_change >= num_iter_zero_change_max ;
		break_isnan_energy = std::isnan( energy ) ;
		
		if( break_change_objective ) {
			break ;
		} else if( break_step_size ) {
			break ;
		} else if( break_sum_net_force_magnitude ) {
			break ;
		} else if( break_zero_change ) {
			break ;
		} else if( break_isnan_energy ) {
			break ;
		}
	}
	std::cout
		<< std::setprecision(3)
		<< std::scientific
		<< "  " << std::setw(10) << energy
		<< "  " << std::setw(10) << energy_prev
		<< "  " << std::setw(10) << energy-energy_prev
		<< "  " << std::setw(10) << step_size
		<< "  " << std::setw(10) << sum_net_force_magnitude
		<< std::endl ;

	std::cout << "solve time: " << omp_get_wtime() - time_start << std::endl ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███    ███  ██████  ██    ██ ███████         ██████   ██████  ██ ███    ██ ████████ ███████         ███    ██ ███████ ██     ██ ████████  ██████  ███    ██ 
// ████  ████ ██    ██ ██    ██ ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ████   ██ ██      ██     ██    ██    ██    ██ ████   ██ 
// ██ ████ ██ ██    ██ ██    ██ █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         ██ ██  ██ █████   ██  █  ██    ██    ██    ██ ██ ██  ██ 
// ██  ██  ██ ██    ██  ██  ██  ██              ██      ██    ██ ██ ██  ██ ██    ██         ██         ██  ██ ██ ██      ██ ███ ██    ██    ██    ██ ██  ██ ██ 
// ██      ██  ██████    ████   ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ██   ████ ███████  ███ ███     ██     ██████  ██   ████ 
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████          ██████  ██████   █████  ██████  ██ ███████ ███    ██ ████████ 
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ██       ██   ██ ██   ██ ██   ██ ██ ██      ████   ██    ██    
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ██   ███ ██████  ███████ ██   ██ ██ █████   ██ ██  ██    ██    
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██    ██ ██   ██ ██   ██ ██   ██ ██ ██      ██  ██ ██    ██    
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████  ██████  ██   ██ ██   ██ ██████  ██ ███████ ██   ████    ██    
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████         ██   ██ ███████ ███████ ███████ ██  █████  ███    ██         ███    ██ ██    ██ ███    ███ ███████ ██████  ██  ██████  █████  ██      
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ██   ██ ██      ██      ██      ██ ██   ██ ████   ██         ████   ██ ██    ██ ████  ████ ██      ██   ██ ██ ██      ██   ██ ██      
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ███████ █████   ███████ ███████ ██ ███████ ██ ██  ██         ██ ██  ██ ██    ██ ██ ████ ██ █████   ██████  ██ ██      ███████ ██      
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██   ██ ██           ██      ██ ██ ██   ██ ██  ██ ██         ██  ██ ██ ██    ██ ██  ██  ██ ██      ██   ██ ██ ██      ██   ██ ██      
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████ ██   ██ ███████ ███████ ███████ ██ ██   ██ ██   ████ ███████ ██   ████  ██████  ██      ██ ███████ ██   ██ ██  ██████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_hessian_numerical( void )
{
	// evaluate spatial derivatives in nodal forces
	// using central difference, small step +/- in each direction
	// select only springs connected to node n to re-evaluate forces
	T delta_position = scale_length * static_cast<T>(1.0e-4) ;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████         ██   ██ ███████ ███████ ███████ ██  █████  ███    ██          █████  ███    ██  █████  ██      ██    ██ ████████ ██  ██████  █████  ██      
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ██   ██ ██      ██      ██      ██ ██   ██ ████   ██         ██   ██ ████   ██ ██   ██ ██       ██  ██     ██    ██ ██      ██   ██ ██      
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ███████ █████   ███████ ███████ ██ ███████ ██ ██  ██         ███████ ██ ██  ██ ███████ ██        ████      ██    ██ ██      ███████ ██      
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██   ██ ██           ██      ██ ██ ██   ██ ██  ██ ██         ██   ██ ██  ██ ██ ██   ██ ██         ██       ██    ██ ██      ██   ██ ██      
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████ ██   ██ ███████ ███████ ███████ ██ ██   ██ ██   ████ ███████ ██   ██ ██   ████ ██   ██ ███████    ██       ██    ██  ██████ ██   ██ ███████ 
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████         ███    ██ ███████ ██     ██ ████████  ██████  ███    ██         ███████ ████████ ███████ ██████          ██████  ██ ██████  ███████  ██████ ████████ ██  ██████  ███    ██ 
// ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ████   ██ ██      ██     ██    ██    ██    ██ ████   ██         ██         ██    ██      ██   ██         ██   ██ ██ ██   ██ ██      ██         ██    ██ ██    ██ ████   ██ 
// ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ██ ██  ██ █████   ██  █  ██    ██    ██    ██ ██ ██  ██         ███████    ██    █████   ██████          ██   ██ ██ ██████  █████   ██         ██    ██ ██    ██ ██ ██  ██ 
// ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██  ██ ██ ██      ██ ███ ██    ██    ██    ██ ██  ██ ██              ██    ██    ██      ██              ██   ██ ██ ██   ██ ██      ██         ██    ██ ██    ██ ██  ██ ██ 
//  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████ ██   ████ ███████  ███ ███     ██     ██████  ██   ████ ███████ ███████    ██    ███████ ██      ███████ ██████  ██ ██   ██ ███████  ██████    ██    ██  ██████  ██   ████ 
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

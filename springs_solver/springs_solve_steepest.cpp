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
// ███    ███ ██ ███    ██ ██ ███    ███ ██ ███████ ███████         ███████ ███    ██ ███████ ██████   ██████  ██    ██ 
// ████  ████ ██ ████   ██ ██ ████  ████ ██    ███  ██              ██      ████   ██ ██      ██   ██ ██        ██  ██  
// ██ ████ ██ ██ ██ ██  ██ ██ ██ ████ ██ ██   ███   █████           █████   ██ ██  ██ █████   ██████  ██   ███   ████   
// ██  ██  ██ ██ ██  ██ ██ ██ ██  ██  ██ ██  ███    ██              ██      ██  ██ ██ ██      ██   ██ ██    ██    ██    
// ██      ██ ██ ██   ████ ██ ██      ██ ██ ███████ ███████ ███████ ███████ ██   ████ ███████ ██   ██  ██████     ██    
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::minimize_energy( void )
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
	std::size_t num_iter_zero_change_max = 100 ;
	T step_size = 1E-3 ;
	T step_size_reduction = 0.9 ; // in range (0,1)
	T step_size_min = 1E-20 ;
	T change_energy ;

	//
	bool is_shaking = false ;
	std::size_t num_iter_shake = 100 ;
	T shake_step_size = 1E-2 ;
	T shake_step_size_reduction = 0.999 ;

	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 500000 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :      0 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :    200 ;
	T local_tolerance_sum_net_force    = ( tolerance_sum_net_force    > 0.0 ) ? tolerance_sum_net_force : 1E-12 ;
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
	//#pragma omp parallel
	for( std::size_t iter = 0 ; iter < local_num_iter_max ; ++iter ) {

		// perturb the system with some randomness
		if( is_shaking ) {
			if( iter%num_iter_shake == 0 ) {
				move_points_rand( shake_step_size , rng , uni_1_1 ) ;
				shake_step_size *= shake_step_size_reduction ;
				update_springs() ;
				update_forces() ; // also computes sum_net_force_magnitude
			}
		}

		// basic line search
		energy_prev = energy ;
		points_prev = points ;
		step_size /= step_size_reduction ;
		step_size /= step_size_reduction ;
		do {
			step_size *= step_size_reduction ;
			points = points_prev ;
			move_points_force( step_size ) ;
			update_springs() ;
			update_forces() ; // also computes sum_net_force_magnitude
			energy = get_objective() ;
		} while( (energy>energy_prev) && (step_size>step_size_min) ) ;

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

		//
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

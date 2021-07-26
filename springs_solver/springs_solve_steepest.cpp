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
// ███    ███ ██ ███    ██ ██ ███    ███ ██ ███████ ███████         ███████ ███    ██ ███████ ██████   ██████  ██    ██ 
// ████  ████ ██ ████   ██ ██ ████  ████ ██    ███  ██              ██      ████   ██ ██      ██   ██ ██        ██  ██  
// ██ ████ ██ ██ ██ ██  ██ ██ ██ ████ ██ ██   ███   █████           █████   ██ ██  ██ █████   ██████  ██   ███   ████   
// ██  ██  ██ ██ ██  ██ ██ ██ ██  ██  ██ ██  ███    ██              ██      ██  ██ ██ ██      ██   ██ ██    ██    ██    
// ██      ██ ██ ██   ████ ██ ██      ██ ██ ███████ ███████ ███████ ███████ ██   ████ ███████ ██   ██  ██████     ██    
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::minimize_energy( void )
{
	// initial state
	update_springs( springs_used ) ;
	update_forces() ;
	energy = total_energy() ;
	get_net_force_mag() ;
	T obj = get_objective() ;
	T obj_prev = obj ;
	T obj_init = obj ;
	//
	std::size_t num_iter_zero_change = 0 ;
	std::size_t num_iter_zero_change_max = 100 ;
	T step_size_reduction = 0.9 ; // in range (0,1)
	T step_size_min = 1E-20 ;
	T change_obj ;
	T mean_obj_change = static_cast<T>(0.0) ;
	T step_size_shared_max ;
	//
	bool is_shaking = false ;
	std::size_t num_iter_shake = 100 ;
	T shake_step_size = 1E-2 ;
	T shake_step_size_reduction = 0.999 ;
	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 5000 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :    0 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :  100 ;
	//
	std::size_t iter_parallel_max = 100 ;
	//
	std::size_t iter = 0 ;
	//
	T local_tolerance_sum_net_force    = ( tolerance_sum_net_force    > 0.0 ) ? tolerance_sum_net_force : 1E-12 ;
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
		std::uniform_real_distribution<T> uni_1_1_local = uni_1_1 ;
		KleinSummer<T> ksum_local ;
		T total_energy_local ;
		T max_net_force_magnitude_local ;
		T sum_net_force_magnitude_local ;
		T step_size_local ;
		// begin gradient descent
		while( !break_any ) {
			// perturb the system with some randomness
			if( is_shaking ) {
				if( iter%num_iter_shake == 0 ) {
					move_points_rand( shake_step_size , rng_local , uni_1_1_local , points_local ) ;
					shake_step_size *= shake_step_size_reduction ;
					update_points_shared( points_update ) ;
					update_springs( springs_local ) ;
					update_forces( points_local , springs_local ) ;
				}
			}
			#pragma omp barrier
			get_net_force_mag( points_local , ksum_local , max_net_force_magnitude_local , sum_net_force_magnitude_local ) ;
			total_energy_local = total_energy( springs_local , thread_springs[thread_id] , ksum_local ) ;
			T obj_local = get_objective( total_energy_local , sum_net_force_magnitude_local , max_net_force_magnitude_local ) ;
			for( std::size_t iter_parallel = 0 ; iter_parallel < iter_parallel_max ; ++iter_parallel ) {
				// basic line search
				T obj_prev_local = obj_local ;
				step_size_local = 1E-3 / (step_size_reduction*step_size_reduction) ;
				points_local_prev = points_local ;
				do {
					step_size_local *= step_size_reduction ;
					points_local = points_local_prev ;
					move_points_force( step_size_local , points_local ) ;
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

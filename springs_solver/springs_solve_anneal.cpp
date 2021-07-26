#include "springs.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
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
//  █████  ███    ██ ███    ██ ███████  █████  ██      
// ██   ██ ████   ██ ████   ██ ██      ██   ██ ██      
// ███████ ██ ██  ██ ██ ██  ██ █████   ███████ ██      
// ██   ██ ██  ██ ██ ██  ██ ██ ██      ██   ██ ██      
// ██   ██ ██   ████ ██   ████ ███████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::anneal( void )
{
	// initial state
	update_springs( springs_used ) ;
	update_forces() ;
	energy = total_energy() ;
	get_net_force_mag() ;
	T obj = get_objective() ;
	T obj_prev = obj ;
	T obj_best = obj ;
	
	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 5000 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :    0 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :   50 ;
	T local_tolerance_sum_net_force    = ( tolerance_sum_net_force    > 0.0 ) ? tolerance_sum_net_force    : 1E-12 ;
	T local_tolerance_change_objective = ( tolerance_change_objective > 0.0 ) ? tolerance_change_objective : 1E-12 ;
	//
	std::size_t iter_parallel_max = 100 ;
	//
	std::size_t iter = 0 ;

	// annealing solver parameters
	// TODO // get these parameters from NetworkParameters? (default/file)
	T step_size = static_cast<T>(0.2) ;
	T step_size_reduction = 0.9 ;
	T temperature_reduction = static_cast<T>(0.99) ;
	T obj_compare ;
	T obj_compare_min = static_cast<T>(1E-12) ;
	T relative_change_obj ;
	T mean_obj_change = static_cast<T>(0.0) ;
	T sum_net_force_magnitude_prev = sum_net_force_magnitude ;
	bool allow_random = false ;
	//std::size_t num_iter_heatup = local_num_iter_max / 5 ;
	std::size_t num_consecutive_reject     =   0 ;
	std::size_t num_consecutive_reject_max = 200 ;
	std::size_t num_since_last_best     = 0 ;
	std::size_t num_since_last_best_max = local_num_iter_max / 100 ;
	std::size_t num_small_change     =  0 ;
	std::size_t num_small_change_max = 20 ;
	std::size_t num_temperature_reductions     =  0 ;
	std::size_t num_temperature_reductions_max = 30 ;
	bool reboot = false ;
	bool reset_to_best = false ;
	bool reset_to_prev = false ;
	bool reset_to_init = false ;
	bool improved_best = false ;
	bool improved_prev = false ;
	bool break_num_iter = false ;
	bool break_temperature_reductions = false ;
	bool break_sum_net_force_magnitude = false ;
	bool break_consecutive_reject = false ;
	bool break_isnan_obj = false ;
	bool break_any = false ;
	find_max_spring_length() ;
	
	// choose starting temperature
	// equivalent to 1/e probability to accept objective difference
	T heatup_amplitude = 0.01 ;
	std::size_t heatup_num_config = 1000 ;
	double time_heat = omp_get_wtime() ;
	obj_compare = heat_up( heatup_amplitude , heatup_num_config ) ;
	std::cout << "heatup time: " << omp_get_wtime() - time_heat << std::endl ;
	obj_compare = ( obj_compare > obj_compare_min ) ? obj_compare : obj_compare_min ;
	T temperature = std::abs( obj_compare - obj_best ) ;

	//
	std::cout
		<< "  " << std::setw(10) << "E_best"
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
	T total_energy_shared           [ num_threads ] ;
	T max_net_force_magnitude_shared[ num_threads ] ;
	T sum_net_force_magnitude_shared[ num_threads ] ;
	bool reboot_shared              [ num_threads ] ;
	#pragma omp parallel
	{
		// assign subsets of modifiable data private to each thread
		const int thread_id = omp_get_thread_num() ;
		std::vector< Point<T,N> > points_local ;
		std::vector< Spring<T,N> > springs_local ;
		std::vector< std::vector< Link > > links_local ;
		std::vector< std::pair< Point<T,N> * , Point<T,N> * > > points_update ;
		setup_local_data( thread_id , points_local , springs_local , links_local , points_update ) ;
		std::vector< Point<T,N> > points_local_init = points_local ;
		std::vector< Point<T,N> > points_local_prev = points_local ;
		std::vector< Point<T,N> > points_local_best = points_local ;
		std::default_random_engine rng_local{ std::random_device{}() } ;
		std::uniform_real_distribution<T> uni_0_1_local = uni_0_1 ;
		std::uniform_real_distribution<T> uni_1_1_local = uni_1_1 ;
		KleinSummer<T> ksum_local ;
		T total_energy_local ;
		T max_net_force_magnitude_local ;
		T sum_net_force_magnitude_local ;
		// begin annealing
		while( !break_any ) {
			//
			     if( reboot        ) { points_local      = points_local_best ; }
			else if( reset_to_init ) { points_local      = points_local_init ; }
			else if( reset_to_prev ) { points_local      = points_local_prev ; }
			else if( reset_to_best ) { points_local      = points_local_best ; points_local_prev = points_local_best ; }
			else if( improved_best ) { points_local_best = points_local      ; points_local_prev = points_local      ; }
			else if( improved_prev ) { points_local_prev = points_local      ; }
			// test a new configuration, force-driven or random
			// apply many sub-iterations of force-driven steps with small step sizes
			// no barriers! asynchronous updates allowed
			if( allow_random ) {
				move_points_rand( step_size * 1E-5 , rng_local , uni_1_1_local , points_local ) ;
			}
			for( std::size_t iter_parallel = 0 ; iter_parallel < iter_parallel_max ; ++iter_parallel ) {
				update_points_shared( points_update ) ;
				update_springs( springs_local ) ;
				update_forces( points_local , springs_local ) ;
				move_points_force( step_size , points_local ) ;
			}
			#pragma omp barrier
			update_points_shared( points_update ) ;
			update_springs( springs_local ) ;
			update_forces( points_local , springs_local ) ;
			get_net_force_mag( points_local , ksum_local , max_net_force_magnitude_local , sum_net_force_magnitude_local ) ;
			total_energy_local = total_energy( springs_local , thread_springs[thread_id] , ksum_local ) ;
			total_energy_shared           [ thread_id ] = total_energy_local            ;
			max_net_force_magnitude_shared[ thread_id ] = max_net_force_magnitude_local ;
			sum_net_force_magnitude_shared[ thread_id ] = sum_net_force_magnitude_local ;
			reboot_shared                 [ thread_id ] = test_reboot( springs_local ) ;
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
				reboot = std::any_of( &reboot_shared[0] ,
									  &reboot_shared[0] + num_threads ,
									  [](const bool & tf){ return tf ; } ) ;
				// evaluate objective function at new configuration
				obj = get_objective() ;
				// if configuration is too extreme, reset network to best
				reset_to_best = false ;
				reset_to_prev = false ;
				reset_to_init = false ;
				improved_best = false ;
				improved_prev = false ;
				if( reboot ) {
					std::cout << "REBOOT" << std::endl ;
					reset_to_best = true ;
					obj_prev = obj_best ;
					reboot = false ;
				} else {

					//
					relative_change_obj = std::abs( obj - obj_prev ) / obj_compare ;
					if( relative_change_obj < local_tolerance_change_objective ) {
						++num_small_change ;
					}

					// accept or reject new state based on change in obj
					// accepting decreased obj is guaranteed
					// accepting increased obj is more likely at high temperatures
					if( obj < obj_best ) {
						improved_best = true ;
						mean_obj_change += (obj-obj_prev) ;
						step_size /= step_size_reduction ;
						obj_best = obj ;
						obj_prev = obj ;
						num_consecutive_reject = 0 ;
						num_since_last_best = 0 ;
					} else if( (obj==obj_best) && (sum_net_force_magnitude<sum_net_force_magnitude_prev) ) {
						// why is the final sum_F_net higher than the one just before it????
						improved_best = true ;
						obj_prev = obj ;
						num_consecutive_reject = 0 ;
						num_since_last_best = 0 ;
					} else if( accept_new_points( obj-obj_prev , temperature , uni_0_1_local(rng_local) ) ) {
						improved_prev = true ;
						mean_obj_change += (obj-obj_prev) ;
						obj_prev = obj ;
						num_consecutive_reject = 0 ;
						++num_since_last_best ;
					} else {
						reset_to_prev = true ;
						step_size *= step_size_reduction ;
						++num_consecutive_reject ;
						++num_since_last_best ;
					}
					//
					if( (local_num_iter_print>0) && (iter%local_num_iter_print==0) ) {
						mean_obj_change /= static_cast<T>(local_num_iter_print) ;
						time_curr = omp_get_wtime() ;
						std::cout
							<< std::setprecision(3)
							<< std::scientific
							<< "  " << std::setw(10) << obj_best
							<< "  " << std::setw(10) << obj
							<< "  " << std::setw(10) << obj_prev
							<< "  " << std::setw(10) << mean_obj_change
							<< "  " << std::setw(10) << step_size
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

					// if obj changes are small for multiple iterations, then
					// assume the current state is near local minimum, and
					// reduce the temperature
					if(    ( num_small_change    >= num_small_change_max    )
						|| ( num_since_last_best >= num_since_last_best_max ) )
					{
						temperature *= temperature_reduction ;
						++num_temperature_reductions ;
						if( num_temperature_reductions >= num_temperature_reductions_max ) {
							break_temperature_reductions = true ;
						}
						num_small_change = 0 ;
						num_since_last_best = 0 ;
					}

					// stopping conditions
					++iter ;
					break_num_iter = iter >= local_num_iter_max ;
					break_sum_net_force_magnitude = sum_net_force_magnitude < local_tolerance_sum_net_force ;
					break_consecutive_reject = num_consecutive_reject >= num_consecutive_reject_max ;
					break_isnan_obj = std::isnan( obj ) ;
					break_any = break_num_iter ||
								break_temperature_reductions ||
								break_sum_net_force_magnitude ||
								break_consecutive_reject ||
								break_isnan_obj ;
				}
			}

			if( reboot ) {
				continue ;
			} else if( break_any ) {
				// synchronize best point positions to shared data
				#pragma omp critical
				{
					if( improved_best ) {
						update_points_all( points_local ) ;
					} else {
						update_points_all( points_local_best ) ;
					}
				}
				break ;
			}

			// FIGURE OUT ANOTHER WAY TO ACHIEVE HEAT UP USING PARALLEL?
			/*
			// if no temperature reductions after many iterations, reset and
			// increase chance of accepting higher obj configurations
			if( (num_temperature_reductions==0) && ( ((iter+1)%num_iter_heatup)==0) ) {
				points = points_init ;
				obj = obj_init ;
				switch( iter / num_iter_heatup ) {
				case 1: heatup_amplitude = 0.01 ; heatup_num_config = 1000 ; break ;
				case 2: heatup_amplitude = 0.01 ; heatup_num_config = 1000 ; break ;
				case 3: heatup_amplitude = 0.02 ; heatup_num_config = 2000 ; break ;
				case 4: heatup_amplitude = 0.05 ; heatup_num_config = 5000 ; break ;
				}
				obj_best = heat_up( heatup_amplitude , heatup_num_config ) ;
				iter += heatup_num_config ;
				points_prev = points ;
				obj_prev = obj ;
			}
			//*/
		}
	}
	std::cout
		<< std::setprecision(3)
		<< std::scientific
		<< "  " << std::setw(10) << obj_best
		<< "  " << std::setw(10) << obj
		<< "  " << std::setw(10) << obj_prev
		<< "  " << std::setw(10) << obj-obj_prev
		<< "  " << std::setw(10) << step_size
		<< "  " << std::setw(10) << sum_net_force_magnitude
		<< std::endl ;

	std::cout << obj - obj_best << std::endl ;

	// final update in best configuration
	update_springs() ;
	update_forces() ;
	energy = total_energy() ;
	get_net_force_mag() ;
	obj = get_objective() ;

	std::cout
		<< std::setprecision(3)
		<< std::scientific
		<< "  " << std::setw(10) << obj_best
		<< "  " << std::setw(10) << obj
		<< "  " << std::setw(10) << obj_prev
		<< "  " << std::setw(10) << obj-obj_prev
		<< "  " << std::setw(10) << step_size
		<< "  " << std::setw(10) << sum_net_force_magnitude
		<< std::endl ;

	std::cout << "solve time: " << omp_get_wtime() - time_start << std::endl ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  █████   ██████  ██████ ███████ ██████  ████████         ███    ██ ███████ ██     ██         ██████   ██████  ██ ███    ██ ████████ ███████ 
// ██   ██ ██      ██      ██      ██   ██    ██            ████   ██ ██      ██     ██         ██   ██ ██    ██ ██ ████   ██    ██    ██      
// ███████ ██      ██      █████   ██████     ██            ██ ██  ██ █████   ██  █  ██         ██████  ██    ██ ██ ██ ██  ██    ██    ███████ 
// ██   ██ ██      ██      ██      ██         ██            ██  ██ ██ ██      ██ ███ ██         ██      ██    ██ ██ ██  ██ ██    ██         ██ 
// ██   ██  ██████  ██████ ███████ ██         ██    ███████ ██   ████ ███████  ███ ███  ███████ ██       ██████  ██ ██   ████    ██    ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
bool SpringNetwork<T,N>::accept_new_points( const T & delta_obj , const T & temperature , const T rand_0_1 )
{
	if( delta_obj < static_cast<T>(0) ) {
		return true ;
	}
	T prob = std::exp( -delta_obj / temperature ) ;
	return ( rand_0_1 < prob ) ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██   ██ ███████  █████  ████████         ██    ██ ██████  
// ██   ██ ██      ██   ██    ██            ██    ██ ██   ██ 
// ███████ █████   ███████    ██            ██    ██ ██████  
// ██   ██ ██      ██   ██    ██            ██    ██ ██      
// ██   ██ ███████ ██   ██    ██    ███████  ██████  ██      
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T SpringNetwork<T,N>::heat_up( const T & amplitude , const std::size_t & num_config_test )
{
	T obj = static_cast<T>(0) ;
	T obj_max = std::numeric_limits<T>::min() ;
	points_prev = points ;
	for( std::size_t i = 0 ; i < num_config_test ; ++i ) {
		move_points_rand( amplitude ) ;
		update_springs() ;
		update_forces() ;
		energy = total_energy() ;
		get_net_force_mag() ;
		obj = get_objective() ;
		if( obj > obj_max ) {
			obj_max = obj ;
		}
		points = points_prev ;
	}
	return obj_max ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ████████ ███████ ███████ ████████         ██████  ███████ ██████   ██████   ██████  ████████ 
//    ██    ██      ██         ██            ██   ██ ██      ██   ██ ██    ██ ██    ██    ██    
//    ██    █████   ███████    ██            ██████  █████   ██████  ██    ██ ██    ██    ██    
//    ██    ██           ██    ██            ██   ██ ██      ██   ██ ██    ██ ██    ██    ██    
//    ██    ███████ ███████    ██    ███████ ██   ██ ███████ ██████   ██████   ██████     ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
bool SpringNetwork<T,N>::test_reboot( void )
{
	iterSpring s ;
	iterSpring s_end ;
	for( iterSpring s = springs.begin() , s_end = springs.end() ; s != s_end ; ++s ) {
		if( s->length > max_spring_length ) {
			return true ;
		}
	}
	return false ;
}

template< class T , std::size_t N >
bool SpringNetwork<T,N>::test_reboot( const std::vector< Spring<T,N> > & springs_subset )
{
	for( std::size_t s = 0 ; s < springs_subset.size() ; ++s ) {
		if( springs_subset[s].length > max_spring_length ) {
			return true ;
		}
	}
	return false ;
}

template< class T , std::size_t N >
bool SpringNetwork<T,N>::test_reboot( const std::vector< Spring<T,N> * > & springs_subset )
{
	for( std::size_t s = 0 ; s < springs_subset.size() ; ++s ) {
		if( springs_subset[s]->length > max_spring_length ) {
			return true ;
		}
	}
	return false ;
}
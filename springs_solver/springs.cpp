//
//  springs.cpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#include "springs.hpp"
#include "vectors_nd.hpp"

#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <random>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>
#include <cstdint>

#include <sys/stat.h>

///______________________ make_dir ______________________///
void make_dir( const std::string & dir_name )
{
	struct stat info ;
	if( stat( dir_name.c_str() , &info ) != 0 ) {
		std::system( ("mkdir " + dir_name).c_str() ) ;
	} else if( info.st_mode & S_IFDIR ) {
		// directory already exists
	} else {
		// exists but not a directory
	}
	return ;
}

////////////////////////////////////////////////////////////
///////////////////// NETWORK PARAMETERS ///////////////////
////////////////////////////////////////////////////////////

NetworkParameters::NetworkParameters( const std::string & dir_input ,
									  const std::string & dir_output )
{
	this->dir_input        = dir_input ;
	this->dir_output       = dir_output ;
	make_dir( dir_output ) ;
	file_input_parameters  = dir_input  + "network_parameters.txt" ;
	file_output_parameters = dir_output + "network_parameters.txt" ;
	load_parameters( file_input_parameters.c_str() ) ;
	save_parameters( file_output_parameters.c_str() ) ;
	return ;
}


void NetworkParameters::load_parameters( const char * file_name )
{
	std::ifstream file ;
	if( file_name != NULL ) {
		file.open( file_name ) ;
		if( !file ) {
			std::cout << "COULD NOT OPEN FILE:\n" << file_name << std::endl ;
		} else {
			file >> num_points
				 >> num_springs
				 >> precision
				 >> num_dimensions
				 >> num_stiffness_tension
				 >> num_stiffness_compression ;
			file.close() ;
		}
	}
	return ;
}

void NetworkParameters::save_parameters( const char * file_name )
{
	std::ofstream file ;
	if( file_name != NULL ) {
		file.open( file_name ) ;
		if( !file ) {
			std::cout << "COULD NOT OPEN FILE:\n" << file_name << std::endl ;
		} else {
			file << num_points                << '\n'
				 << num_springs               << '\n'
				 << precision                 << '\n'
				 << num_dimensions            << '\n'
				 << num_stiffness_tension     << '\n'
				 << num_stiffness_compression ;
			file.close() ;
		}
	}
	return ;
}


////////////////////////////////////////////////////////////
//////////////////////////  POINT //////////////////////////
////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////
////////////////////////// SPRING //////////////////////////
////////////////////////////////////////////////////////////

///___________________  spring_energy ___________________///
template< class T , std::size_t N >
T Spring<T,N>::spring_energy( void )
{
	Vector<T,N> delta_position( end->position ) ;
	delta_position -= start->position ;
	length = delta_position.norm() ;
	T delta_length = length - rest_length ;
	T force_magnitude = static_cast<T>(0) ;
	T energy = static_cast<T>(0) ;
	if( delta_length > 0 ) {
		T force_component ;
		T delta_length_power_i = delta_length ;
		for( std::size_t i = 0 ; i < Spring<T,N>::num_stiffness_tension ; ++i ) {
			force_component = stiffness_tension[i] * delta_length_power_i ;
			force_magnitude += force_component ;
			energy += force_component * delta_length / static_cast<T>(i+1) ;
			delta_length_power_i *= delta_length ;
		}
	} else if( allow_compression ) {
		delta_length = std::fabs( delta_length ) ;
		T force_component ;
		T delta_length_power_i = delta_length ;
		for( std::size_t i = 0 ; i < Spring<T,N>::num_stiffness_compression ; ++i ) {
			force_component = stiffness_compression[i] * delta_length_power_i ;
			force_magnitude += force_component ;
			energy += force_component * delta_length / static_cast<T>(i+1) ;
			delta_length_power_i *= delta_length ;
		}
		force_magnitude = -force_magnitude ;
	}
	force = delta_position ;
	force *= (force_magnitude/length) ;
	return energy ;
}

/// DECLARE STATIC MEMBER VARIABLES ///
template< class T , std::size_t N >
std::size_t Spring<T,N>::num_stiffness_tension ;
template< class T , std::size_t N >
std::size_t Spring<T,N>::num_stiffness_compression ;

////////////////////////////////////////////////////////////
////////////////////// SPRING NETWORK //////////////////////
////////////////////////////////////////////////////////////

///_____________  create_spring_network_obj _____________///
std::unique_ptr<ASpringNetwork> ASpringNetwork::create_spring_network_obj( const NetworkParameters & network_parameters )
{
	std::unique_ptr<ASpringNetwork> asn = NULL ;
	if( network_parameters.precision.compare( "double" ) == 0 ) {
		if     ( network_parameters.num_dimensions == 2 ) { asn = std::unique_ptr<ASpringNetwork>( new SpringNetwork<double,2>() ) ; }
		else if( network_parameters.num_dimensions == 3 ) { asn = std::unique_ptr<ASpringNetwork>( new SpringNetwork<double,3>() ) ; }
		else if( network_parameters.num_dimensions == 4 ) { asn = std::unique_ptr<ASpringNetwork>( new SpringNetwork<double,4>() ) ; }
		else { return asn ; }
	} else if( network_parameters.precision.compare( "float" ) == 0 ) {
		if     ( network_parameters.num_dimensions == 2 ) { asn = std::unique_ptr<ASpringNetwork>( new SpringNetwork<float,2>() ) ; }
		else if( network_parameters.num_dimensions == 3 ) { asn = std::unique_ptr<ASpringNetwork>( new SpringNetwork<float,3>() ) ; }
		else if( network_parameters.num_dimensions == 4 ) { asn = std::unique_ptr<ASpringNetwork>( new SpringNetwork<float,4>() ) ; }
		else { return asn ; }
	} else {
		return asn ;
	}
	return asn ;
}

///_______________________  setup _______________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::setup( const NetworkParameters & network_parameters )
{
	dir_input           = network_parameters.dir_input ;
	dir_output          = network_parameters.dir_output ;
	dir_output_iter     = dir_output + "iter/" ;
	if( num_iter_save > 0 ) {
		make_dir( dir_output_iter ) ;
	}
	file_input_nodes    = dir_input  + "network_nodes.dat" ;
	file_input_springs  = dir_input  + "network_springs.dat" ;
	file_output_nodes   = dir_output + "network_nodes.dat" ;
	file_output_springs = dir_output + "network_springs.dat" ;
	uni_0_1 = std::uniform_real_distribution<T>( 0.0,+1.0) ;
	uni_1_1 = std::uniform_real_distribution<T>(-1.0,+1.0) ;
	num_points  = network_parameters.num_points ;
	num_springs = network_parameters.num_springs ;
	Spring<T,N>::num_stiffness_tension = network_parameters.num_stiffness_tension ;
	Spring<T,N>::num_stiffness_compression = network_parameters.num_stiffness_compression ;
	points  = std::vector<  Point<T,N> >( num_points  ) ;
	springs = std::vector< Spring<T,N> >( num_springs ) ;
	nodes   = std::vector<        Node >( num_points  ) ;
	load_network_binary( file_input_nodes.c_str() , file_input_springs.c_str() ) ;
	construct_network() ;
	return ;
}

///____________________  read_binary ____________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::load_network_binary( const char * file_nodes ,
											  const char * file_springs )
{
	std::FILE * file_ptr ;
	T             * read_data_T ;
	std::uint32_t * read_data_uint32 ;
	std::uint8_t  * read_data_uint8 ;
	
	std::size_t NKT = Spring<T,N>::num_stiffness_tension ;
	std::size_t NKC = Spring<T,N>::num_stiffness_compression ;
	
	// points
	if( file_nodes != NULL ) {
		file_ptr = std::fopen( file_nodes , "rb" ) ;
		read_data_T     = new T            [N+N] ;
		read_data_uint8 = new std::uint8_t [1  ] ;
		for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
			std::fread( read_data_T     , sizeof(T)            , N+N , file_ptr ) ;
			std::fread( read_data_uint8 , sizeof(std::uint8_t) , 1   , file_ptr ) ;
			for( std::size_t n = 0 ; n < N ; ++n ) {
				p->position[n]      = read_data_T[n  ] ;
				p->force_applied[n] = read_data_T[n+N] ;
			}
			p->fixed = ( read_data_uint8[0] == 1 ) ;
		}
		delete [] read_data_T     ;
		delete [] read_data_uint8 ;
		std::fclose( file_ptr ) ;
	}
	
	// springs
	if( file_springs != NULL ) {
		file_ptr = std::fopen( file_springs , "rb" ) ;
		read_data_uint32 = new std::uint32_t [2        ] ;
		read_data_T      = new T             [1+NKT+NKC] ;
		read_data_uint8  = new std::uint8_t  [1        ] ;
		for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
			std::fread( read_data_uint32 , sizeof(std::uint32_t) , 2         , file_ptr ) ;
			std::fread( read_data_T      , sizeof(T)             , 1+NKT+NKC , file_ptr ) ;
			std::fread( read_data_uint8  , sizeof(std::uint8_t)  , 1         , file_ptr ) ;
			s->start = &points[ read_data_uint32[0] ] ;
			s->end   = &points[ read_data_uint32[1] ] ;
			s->rest_length = read_data_T[0] ;
			s->stiffness_tension     = std::vector<T>( &read_data_T[1    ] , &read_data_T[1+NKT    ] ) ;
			s->stiffness_compression = std::vector<T>( &read_data_T[1+NKT] , &read_data_T[1+NKT+NKC] ) ;
			s->allow_compression = ( read_data_uint8[0] == 1 ) ;
		}
		delete [] read_data_uint32 ;
		delete [] read_data_T      ;
		delete [] read_data_uint8  ;
		std::fclose( file_ptr ) ;
	}
	return ;
}

///____________________ write_binary ____________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::save_network_binary( const char * file_nodes ,
											  const char * file_springs )
{
	std::FILE * file_ptr ;
	T             * write_data_T ;
	std::uint32_t * write_data_uint32 ;
	std::uint8_t  * write_data_uint8 ;
	
	std::size_t NKT = Spring<T,N>::num_stiffness_tension ;
	std::size_t NKC = Spring<T,N>::num_stiffness_compression ;
	
	// points
	if( file_nodes != NULL ) {
		file_ptr = std::fopen( file_nodes , "wb" ) ;
		write_data_T     = new T            [N+N] ;
		write_data_uint8 = new std::uint8_t [1  ] ;
		for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
			for( std::size_t n = 0 ; n < N ; ++n ) {
				write_data_T[n  ] = p->position[n]      ;
				write_data_T[n+N] = p->force_applied[n] ;
			}
			write_data_uint8[0] = static_cast<bool>( p->fixed ) ;
			std::fwrite( write_data_T     , sizeof(T)            , N+N , file_ptr ) ;
			std::fwrite( write_data_uint8 , sizeof(std::uint8_t) , 1   , file_ptr ) ;
		}
		delete [] write_data_T     ;
		delete [] write_data_uint8 ;
		std::fclose( file_ptr ) ;
	}
	
	// springs
	if( file_springs != NULL ) {
		file_ptr = std::fopen( file_springs , "wb" ) ;
		write_data_uint32 = new std::uint32_t [2        ] ;
		write_data_T      = new T             [1+NKT+NKC] ;
		write_data_uint8  = new std::uint8_t  [1        ] ;
		for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
			write_data_uint32[0] = static_cast<std::uint32_t>( s->start - &points[0] ) ;
			write_data_uint32[1] = static_cast<std::uint32_t>( s->end   - &points[0] ) ;
			write_data_T[0] = s->rest_length ;
			std::copy( s->stiffness_tension.begin()     , s->stiffness_tension.end()     , &write_data_T[1    ] ) ;
			std::copy( s->stiffness_compression.begin() , s->stiffness_compression.end() , &write_data_T[1+NKT] ) ;
			write_data_uint8[0] = static_cast<bool>( s->allow_compression ) ;
			std::fwrite( write_data_uint32 , sizeof(std::uint32_t) , 2         , file_ptr ) ;
			std::fwrite( write_data_T      , sizeof(T)             , 1+NKT+NKC , file_ptr ) ;
			std::fwrite( write_data_uint8  , sizeof(std::uint8_t)  , 1         , file_ptr ) ;
		}
		delete [] write_data_uint32 ;
		delete [] write_data_T      ;
		delete [] write_data_uint8  ;
		std::fclose( file_ptr ) ;
	}
	return ;
}

///_________________  construct_network _________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::construct_network( void )
{
	for( std::size_t p = 0 ; p < num_points ; ++p ) {
		nodes[p].point = &points[p] ;
	}
	std::size_t p_start ;
	std::size_t p_end ;
	for( std::size_t s = 0 ; s < num_springs ; ++s ) {
		p_start = springs[s].start - &points[0] ;
		p_end   = springs[s].end   - &points[0] ;
		nodes[p_start].links.push_back( (Link){ &points[p_end  ] , &springs[s] , static_cast<T>(+1) } ) ;
		nodes[p_end  ].links.push_back( (Link){ &points[p_start] , &springs[s] , static_cast<T>(-1) } ) ;
	}
	nodes.erase( std::remove_if( nodes.begin() , nodes.end() , []( Node n ){ return n.point->fixed ; } ) , nodes.end() ) ;
	return ;
}

///_______________ find_max_spring_length _______________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::find_max_spring_length( void )
{
	Vector<T,N> x_min ;
	Vector<T,N> x_max ;
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		for( std::size_t n = 0 ; n < N ; ++n ) {
			x_min[n] = ( p->position[n] < x_min[n] ) ? p->position[n] : x_min[n] ;
			x_max[n] = ( p->position[n] > x_max[n] ) ? p->position[n] : x_max[n] ;
		}
	}
	max_spring_length = static_cast<T>(0) ;
	T dx ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		dx = x_max[n] - x_min[n] ;
		max_spring_length = ( dx > max_spring_length ) ? dx : max_spring_length ;
	}
	return ;
}

///____________________ total_energy ____________________///
template< class T , std::size_t N >
T SpringNetwork<T,N>::total_energy( void )
{
	// internal energy due to spring tension/compression
	T energy = static_cast<T>(0) ;
	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
		energy += s->spring_energy() ;
	}
	
	// external work done by applied loads
	Vector<T,N> displacement ;
	for( std::size_t p = 0 ; p < num_points ; ++p ) {
		displacement = points[p].position ;
		displacement -= points_init[p].position ;
		energy -= points[p].force_applied.dot( displacement ) ;
	}
	
	// account for force imbalances, "potential energy"
	// combine external applied force & internal spring forces
	iterNode n ;
	iterNode n_end ;
	iterLink l ;
	iterLink l_end ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		n->net_force = n->point->force_applied ;
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			if( l->spring_direction > static_cast<T>(0) ) {
				n->net_force += l->spring->get_force() ;
			} else {
				n->net_force -= l->spring->get_force() ;
			}
		}
		n->net_force_magnitude = n->net_force.norm() ;
	}
	T sum_spring_length = static_cast<T>(0) ;
	T sum_net_force_magnitude = static_cast<T>(0) ;
	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
		sum_spring_length += s->length ;
	}
	for( iterNode n = nodes.begin() ; n != nodes.end() ; ++n ) {
		sum_net_force_magnitude += n->net_force_magnitude ;
	}
	energy += ( sum_net_force_magnitude * sum_spring_length ) ;
	
	return energy ;
}

///_________________  move_points_force _________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_force( const T & force_step_size )
{
	iterNode n ;
	iterNode n_end ;
	// apply net force & move small displacement towards equilibrating position
	Vector<T,N> small_displacement ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		small_displacement = n->net_force ;
		small_displacement *= force_step_size ;
		n->point->position += small_displacement ;
	}
	return ;
}

///__________________ move_points_rand __________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_rand( const T & amplitude )
{
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		for( std::size_t n = 0 ; n < N ; ++n ) {
			p->position[n] += amplitude * uni_1_1(rng) ;
		}
	}
	return ;
}

///____________________  test_reboot ____________________///
template< class T , std::size_t N >
bool SpringNetwork<T,N>::test_reboot( void )
{
	iterSpring s ;
	iterSpring s_end ;
	for( s = springs.begin() , s_end = springs.end() ; s != s_end ; ++s ) {
		if( s->length > max_spring_length ) {
			return true ;
		}
	}
	return false ;
}

///_________________  accept_new_points _________________///
template< class T , std::size_t N >
bool SpringNetwork<T,N>::accept_new_points( const T & delta_energy , const T & temperature )
{
	if( delta_energy < static_cast<T>(0) ) {
		return true ;
	}
	T rand = uni_0_1(rng) ;
	T prob = std::exp( -delta_energy / temperature ) ;
	return ( rand < prob ) ;
}

///______________________  heat_up ______________________///
template< class T , std::size_t N >
T SpringNetwork<T,N>::heat_up( const T & amplitude , const T & num_config_test )
{
	T energy = static_cast<T>(0) ;
	T energy_max = std::numeric_limits<T>::min() ;
	points_prev = points ;
	for( std::size_t i = 0 ; i < num_config_test ; ++i ) {
		move_points_rand( amplitude ) ;
		energy = total_energy() ;
		if( energy > energy_max ) {
			energy_max = energy ;
		}
		points = points_prev ;
	}
	return energy_max ;
}

///_______________________  solve _______________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::solve( void )
{
	anneal() ;
	save_output() ;
	return ;
}

///____________________  save_output ____________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::save_output( void )
{
	save_network_binary( file_output_nodes.c_str() , file_output_springs.c_str() ) ;
	return ;
}

///_______________________ anneal _______________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::anneal( void )
{
	// copy current state to initial, previous, and best states
	points_init = points ;
	points_prev = points ;
	points_best = points ;
	T energy = total_energy() ;
	T energy_init = energy ;
	T energy_prev = energy ;
	T energy_best = energy ;

	// count how many intermediate iterations saved
	std::size_t num_iter_saved = 0 ;
	
	// annealing solver parameters
	// TODO // get these parameters from NetworkParameters (default/file)
	T force_step_size = static_cast<T>(0.01) ;
	T temperature_reduction = static_cast<T>(0.99) ;
	T energy_compare ;
	T energy_compare_min = static_cast<T>(1E-12) ;
	T relative_change_energy ;
	T relative_change_energy_tol = static_cast<T>(1E-6) ;
	std::size_t num_iter_max       = 500000 ;
	std::size_t num_iter_heatup    = num_iter_max /   5 ;
	std::size_t num_iter_01percent = num_iter_max / 100 ;
	std::size_t num_iter_10percent = num_iter_max /  10 ;
	std::size_t num_small_change     =  0 ;
	std::size_t num_small_change_max = 20 ;
	std::size_t num_temperature_reductions     =    0 ;
	std::size_t num_temperature_reductions_max = 1000 ;
	bool reboot ;
	find_max_spring_length() ;
	
	// choose starting temperature
	// equivalent to 1/e probability to accept energy difference
	T heatup_amplitude = 0.01 ;
	std::size_t heatup_num_config = 1000 ;
	energy_compare = heat_up( heatup_amplitude , heatup_num_config ) ;
	energy_compare = ( energy_compare > energy_compare_min ) ? energy_compare : energy_compare_min ;
	T temperature = std::fabs( energy_compare - energy_best ) ;
	
	// begin annealing
	for( std::size_t iter = 0 ; iter < num_iter_max ; ++iter ) {
		if( iter%num_iter_01percent == 0 ) {
			if( iter > 0 ) {
				std::cout << "." << std::flush ;
				if( iter%num_iter_10percent == 0 ) {
					std::cout << 100*iter/num_iter_max << "%\n\t" << std::flush ;
				}
			} else {
				std::cout << "\t" << std::flush ;
			}
		}
		
		// randomize the order of nodes and test a new configuration
		// if configuration is too extreme, reset network to best
		// no need for random ordering if move_points_force() does not depend on order
		///std::random_shuffle( nodes.begin() , nodes.end() ) ;
		move_points_force( force_step_size ) ;
		energy = total_energy() ;
		reboot = test_reboot() ;
		if( reboot ) {
			std::cout << "REBOOT" << std::endl ;
			points = points_best ;
			points_prev = points_best ;
			energy_prev = energy_best ;
			reboot = false ;
			continue ;
		}
		
		// if energy changes are small for multiple iterations, then
		// assume the current state is near local minimum, and
		// reduce the temperature
		relative_change_energy = std::fabs( energy - energy_prev ) / energy_compare ;
		/**
		if( iter%100 == 0 ) {
			std::cout
				<< "\t" << energy
				<< "\t" << energy_prev
				<< "\t" << relative_change_energy
				<< std::endl ;
		}
		//*/
		if( relative_change_energy < relative_change_energy_tol ) {
			++num_small_change ;
			if( num_small_change == num_small_change_max ) {
				temperature *= temperature_reduction ;
				++num_temperature_reductions ;
				if( num_temperature_reductions == num_temperature_reductions_max ) {
					break ;
				}
				num_small_change = 0 ;
			}
		}
		
		if( num_iter_saved < num_iter_save ) {
			if( ((iter%100)==0) && (iter<=10000) ) {
				std::string dir_output_iter_curr = dir_output_iter + "iter_" + std::to_string(iter) + "/" ;
				make_dir( dir_output_iter_curr ) ;
				save_network_binary(
					(dir_output_iter_curr+"network_nodes.dat").c_str()   ,
					(dir_output_iter_curr+"network_springs.dat").c_str() ) ;
				++num_iter_saved ;
			}
		}
		
		// accept or reject new state based on change in energy
		// accepting decreased energy is guaranteed
		// accepting increased energy is more likely at high temperatures
		if( energy < energy_best ) {
			points_best = points ;
			points_prev = points ;
			energy_best = energy ;
			energy_prev = energy ;
		}
		else if( accept_new_points( energy-energy_prev , temperature ) ) {
			points_prev = points ;
			energy_prev = energy ;
		} else {
			points = points_prev ;
		}
		
		// if no temperature reductions after many iterations, reset and
		// increase chance of accepting higher energy configurations
		if( (num_temperature_reductions==0) && ( ((iter+1)%num_iter_heatup)==0) ) {
			points = points_init ;
			energy = energy_init ;
			switch( iter / num_iter_heatup ) {
			case 1: heatup_amplitude = 0.01 ; heatup_num_config = 1000 ; break ;
			case 2: heatup_amplitude = 0.01 ; heatup_num_config = 1000 ; break ;
			case 3: heatup_amplitude = 0.02 ; heatup_num_config = 2000 ; break ;
			case 4: heatup_amplitude = 0.05 ; heatup_num_config = 5000 ; break ;
			}
			energy_best = heat_up( heatup_amplitude , heatup_num_config ) ;
			iter += heatup_num_config ;
			if( ((iter-heatup_num_config)%num_iter_01percent) >= (num_iter_01percent-heatup_num_config) ) {
				if( iter > 0 ) {
					std::cout << "." << std::flush ;
					if( ((iter-heatup_num_config)%num_iter_10percent) >= (num_iter_10percent-heatup_num_config) ) {
						std::cout << 100*iter/num_iter_max << "%\n\t" << std::flush ;
					}
				} else {
					std::cout << "\t" << std::flush ;
				}
			}
			points_prev = points ;
			energy_prev = energy ;
		}
	}
	points = points_best ;
	energy = energy_best ;
	return ;
}

/// EXPLICIT TEMPLATE INSTANTIATIONS ///
template class Point<float ,2> ;
template class Point<float ,3> ;
template class Point<float ,4> ;
template class Point<double,2> ;
template class Point<double,3> ;
template class Point<double,4> ;

template class Spring<float ,2> ;
template class Spring<float ,3> ;
template class Spring<float ,4> ;
template class Spring<double,2> ;
template class Spring<double,3> ;
template class Spring<double,4> ;

template class SpringNetwork<float ,2> ;
template class SpringNetwork<float ,3> ;
template class SpringNetwork<float ,4> ;
template class SpringNetwork<double,2> ;
template class SpringNetwork<double,3> ;
template class SpringNetwork<double,4> ;



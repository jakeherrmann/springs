//
//  springs.cpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//
//  headers from: http://patorjk.com/software/taag
//  class names "univers"
//  function names "ANSI regular"
//

#include "springs.hpp"
#include "vectors_nd.hpp"
#include "spmat.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <string>

#include <omp.h>
#include <stdlib.h>

#include <sys/stat.h>

#ifdef _WIN32
	#define FILESEP '\\'
#else
	#define FILESEP '/'
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███    ███  █████  ██   ██ ███████         ██████  ██ ██████  
//  ████  ████ ██   ██ ██  ██  ██              ██   ██ ██ ██   ██ 
//  ██ ████ ██ ███████ █████   █████           ██   ██ ██ ██████  
//  ██  ██  ██ ██   ██ ██  ██  ██              ██   ██ ██ ██   ██ 
//  ██      ██ ██   ██ ██   ██ ███████ ███████ ██████  ██ ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  888b      88                                                                    88         88888888ba                                                                                                                   
//  8888b     88                ,d                                                  88         88      "8b                                                                        ,d                                        
//  88 `8b    88                88                                                  88         88      ,8P                                                                        88                                        
//  88  `8b   88   ,adPPYba,  MM88MMM  8b      db      d8   ,adPPYba,   8b,dPPYba,  88   ,d8   88aaaaaa8P'  ,adPPYYba,  8b,dPPYba,  ,adPPYYba,  88,dPYba,,adPYba,    ,adPPYba,  MM88MMM   ,adPPYba,  8b,dPPYba,  ,adPPYba,  
//  88   `8b  88  a8P_____88    88     `8b    d88b    d8'  a8"     "8a  88P'   "Y8  88 ,a8"    88""""""'    ""     `Y8  88P'   "Y8  ""     `Y8  88P'   "88"    "8a  a8P_____88    88     a8P_____88  88P'   "Y8  I8[    ""  
//  88    `8b 88  8PP"""""""    88      `8b  d8'`8b  d8'   8b       d8  88          8888[      88           ,adPPPPP88  88          ,adPPPPP88  88      88      88  8PP"""""""    88     8PP"""""""  88           `"Y8ba,   
//  88     `8888  "8b,   ,aa    88,      `8bd8'  `8bd8'    "8a,   ,a8"  88          88`"Yba,   88           88,    ,88  88          88,    ,88  88      88      88  "8b,   ,aa    88,    "8b,   ,aa  88          aa    ]8I  
//  88      `888   `"Ybbd8"'    "Y888      YP      YP       `"YbbdP"'   88          88   `Y8a  88           `"8bbdP"Y8  88          `"8bbdP"Y8  88      88      88   `"Ybbd8"'    "Y888   `"Ybbd8"'  88          `"YbbdP"'  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██ ██████   █████  ██████   █████  ███    ███ ███████ ████████ ███████ ██████  ███████ 
//  ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██  ██   ██ ██   ██ ██   ██ ██   ██ ████  ████ ██         ██    ██      ██   ██ ██      
//  ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████   ██████  ███████ ██████  ███████ ██ ████ ██ █████      ██    █████   ██████  ███████ 
//  ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██  ██      ██   ██ ██   ██ ██   ██ ██  ██  ██ ██         ██    ██      ██   ██      ██ 
//  ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ ██      ██   ██ ██   ██ ██   ██ ██      ██ ███████    ██    ███████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██       ██████   █████  ██████          ██████   █████  ██████   █████  ███    ███ ███████ ████████ ███████ ██████  ███████ 
//  ██      ██    ██ ██   ██ ██   ██         ██   ██ ██   ██ ██   ██ ██   ██ ████  ████ ██         ██    ██      ██   ██ ██      
//  ██      ██    ██ ███████ ██   ██         ██████  ███████ ██████  ███████ ██ ████ ██ █████      ██    █████   ██████  ███████ 
//  ██      ██    ██ ██   ██ ██   ██         ██      ██   ██ ██   ██ ██   ██ ██  ██  ██ ██         ██    ██      ██   ██      ██ 
//  ███████  ██████  ██   ██ ██████  ███████ ██      ██   ██ ██   ██ ██   ██ ██      ██ ███████    ██    ███████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NetworkParameters::load_parameters( const char * file_name )
{
	std::ifstream file ;
	std::string arg ;
	std::string val ;
	if( file_name != NULL ) {
		file.open( file_name ) ;
		if( !file ) {
			std::cout << "COULD NOT OPEN FILE:\n" << file_name << std::endl ;
		} else {
			while( file >> arg >> val ) {
				if( arg=="num_points"                 ) { num_points                 = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_springs"                ) { num_springs                = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="precision"                  ) { precision                  =                                    val   ; } else
				if( arg=="num_dimensions"             ) { num_dimensions             = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="algorithm"                  ) { algorithm                  =                                    val   ; } else
				if( arg=="objective"                  ) { objective                  =                                    val   ; } else
				if( arg=="num_threads"                ) { user_num_threads           = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_save"              ) { num_iter_save              = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_print"             ) { num_iter_print             = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_max"               ) { num_iter_max               = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="include_force_fixed_nodes"  ) { include_force_fixed_nodes  =        static_cast<bool>(std::stod(val)) ; } else
				if( arg=="use_numerical_hessian"      ) { use_numerical_hessian      =        static_cast<bool>(std::stod(val)) ; } else
				if( arg=="tolerance_change_objective" ) { tolerance_change_objective =                          std::stod(val)  ; } else
				if( arg=="tolerance_sum_net_force"    ) { tolerance_sum_net_force    =                          std::stod(val)  ; }
				else { std::cout << "Parameter not recognized: " << arg << std::endl ; }
			}
			file.close() ;
		}
	}
	return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████  █████  ██    ██ ███████         ██████   █████  ██████   █████  ███    ███ ███████ ████████ ███████ ██████  ███████ 
//  ██      ██   ██ ██    ██ ██              ██   ██ ██   ██ ██   ██ ██   ██ ████  ████ ██         ██    ██      ██   ██ ██      
//  ███████ ███████ ██    ██ █████           ██████  ███████ ██████  ███████ ██ ████ ██ █████      ██    █████   ██████  ███████ 
//       ██ ██   ██  ██  ██  ██              ██      ██   ██ ██   ██ ██   ██ ██  ██  ██ ██         ██    ██      ██   ██      ██ 
//  ███████ ██   ██   ████   ███████ ███████ ██      ██   ██ ██   ██ ██   ██ ██      ██ ███████    ██    ███████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NetworkParameters::save_parameters( const char * file_name )
{
	std::ofstream file ;
	if( file_name != NULL ) {
		file.open( file_name ) ;
		if( !file ) {
			std::cout << "COULD NOT OPEN FILE:\n" << file_name << std::endl ;
		} else {
			file << "num_points"                 << ' ' << num_points                 << '\n'
				 << "num_springs"                << ' ' << num_springs                << '\n'
				 << "precision"                  << ' ' << precision                  << '\n'
				 << "num_dimensions"             << ' ' << num_dimensions             << '\n'
				 << "algorithm"                  << ' ' << algorithm                  << '\n'
				 << "objective"                  << ' ' << objective                  << '\n'
				 << "num_threads"                << ' ' << user_num_threads           << '\n'
				 << "num_iter_save"              << ' ' << num_iter_save              << '\n'
				 << "num_iter_print"             << ' ' << num_iter_print             << '\n'
				 << "num_iter_max"               << ' ' << num_iter_max               << '\n'
				 << "include_force_fixed_nodes"  << ' ' << include_force_fixed_nodes  << '\n'
				 << "use_numerical_hessian"      << ' ' << use_numerical_hessian      << '\n'
				 << "tolerance_change_objective" << ' ' << tolerance_change_objective << '\n'
				 << "tolerance_sum_net_force"    << ' ' << tolerance_sum_net_force ;
			file.close() ;
		}
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  88888888ba                88                        
//  88      "8b               ""                 ,d     
//  88      ,8P                                  88     
//  88aaaaaa8P'   ,adPPYba,   88  8b,dPPYba,   MM88MMM  
//  88""""""'    a8"     "8a  88  88P'   `"8a    88     
//  88           8b       d8  88  88       88    88     
//  88           "8a,   ,a8"  88  88       88    88,    
//  88            `"YbbdP"'   88  88       88    "Y888  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███    ███  ██████  ██    ██ ███████ 
//  ████  ████ ██    ██ ██    ██ ██      
//  ██ ████ ██ ██    ██ ██    ██ █████   
//  ██  ██  ██ ██    ██  ██  ██  ██      
//  ██      ██  ██████    ████   ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Point<T,N>::move( const Vector<T,N> & displacement )
{
	if( !this->fixed_all_dim ) {
		this->position += displacement.zero_where( this->fixed_dim ) ;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   ad88888ba                            88                            
//  d8"     "8b                           ""                            
//  Y8,                                                                 
//  `Y8aaaaa,    8b,dPPYba,   8b,dPPYba,  88  8b,dPPYba,    ,adPPYb,d8  
//    `"""""8b,  88P'    "8a  88P'   "Y8  88  88P'   `"8a  a8"    `Y88  
//          `8b  88       d8  88          88  88       88  8b       88  
//  Y8a     a8P  88b,   ,a8"  88          88  88       88  "8a,   ,d88  
//   "Y88888P"   88`YbbdP"'   88          88  88       88   `"YbbdP"Y8  
//               88                                         aa,    ,88  
//               88                                          "Y8bbdP"   
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████ ███    ██ ███████ ██████   ██████  ██    ██ 
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██      ████   ██ ██      ██   ██ ██        ██  ██  
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         █████   ██ ██  ██ █████   ██████  ██   ███   ████   
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██      ██  ██ ██ ██      ██   ██ ██    ██    ██    
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ███████ ██   ████ ███████ ██   ██  ██████     ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T Spring<T,N>::spring_energy( void )
{
	delta_position = end->position ;
	delta_position -= start->position ;
	length = delta_position.norm() ;
	T delta_length = length - rest_length ;
	T force_magnitude = static_cast<T>(0) ;
	T force_sign = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	Spring<T,N>::ForceLengthRelationship * force_length_type = NULL ;
	std::size_t * num_force_length_parameters = NULL ;
	std::vector<T> * force_length_parameters = NULL ;
	if( delta_length == 0.0 ) {
		force_magnitude = static_cast<T>(0) ;
		energy = static_cast<T>(0) ;
		effective_spring_constant = spring_stiffness_rest() ;
	} else {
		if( delta_length > 0 ) {
			force_length_type           = & force_length_type_tension ;
			num_force_length_parameters = & num_force_length_parameters_tension ;
			force_length_parameters     = & force_length_parameters_tension ;
			force_sign = static_cast<T>(+1) ;
		} else if( delta_length < 0 ) {
			delta_length = std::abs( delta_length ) ;
			force_length_type           = & force_length_type_compression ;
			num_force_length_parameters = & num_force_length_parameters_compression ;
			force_length_parameters     = & force_length_parameters_compression ;
			force_sign = static_cast<T>(-1) ;
		}
		switch( * force_length_type ) {
			case Spring<T,N>::ForceLengthRelationship::polynomial :
				Spring<T,N>::spring_force_polynomial( delta_length ,
													  * num_force_length_parameters ,
													  * force_length_parameters ,
													  effective_spring_constant ,
													  force_magnitude ,
													  energy ) ;
				break ;
			case Spring<T,N>::ForceLengthRelationship::exponential :
				Spring<T,N>::spring_force_exponential( delta_length ,
													   * num_force_length_parameters ,
													   * force_length_parameters ,
													   effective_spring_constant ,
													   force_magnitude ,
													   energy ) ;
				break ;
			case Spring<T,N>::ForceLengthRelationship::powerlaw :
				Spring<T,N>::spring_force_powerlaw( delta_length ,
												    * num_force_length_parameters ,
												    * force_length_parameters ,
													effective_spring_constant ,
												    force_magnitude ,
												    energy ) ;
				break ;
			default:
				effective_spring_constant = static_cast<T>(0) ;
				force_magnitude = static_cast<T>(0) ;
				energy = static_cast<T>(0) ;
		}
		force_magnitude *= force_sign ;
	}
	force = delta_position ;
	force *= (force_magnitude/length) ;
	return energy ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████  ██████  ██████   ██████ ███████         ██████   ██████  ██      ██    ██ ███    ██  ██████  ███    ███ ██  █████  ██      
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██      ██    ██ ██   ██ ██      ██              ██   ██ ██    ██ ██       ██  ██  ████   ██ ██    ██ ████  ████ ██ ██   ██ ██      
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         █████   ██    ██ ██████  ██      █████           ██████  ██    ██ ██        ████   ██ ██  ██ ██    ██ ██ ████ ██ ██ ███████ ██      
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██      ██    ██ ██   ██ ██      ██              ██      ██    ██ ██         ██    ██  ██ ██ ██    ██ ██  ██  ██ ██ ██   ██ ██      
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ██       ██████  ██   ██  ██████ ███████ ███████ ██       ██████  ███████    ██    ██   ████  ██████  ██      ██ ██ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::spring_force_polynomial( const T & delta_length ,
										   const std::size_t & num_force_length_parameters ,
										   const std::vector<T> & force_length_parameters ,
										   T & effective_spring_constant ,
										   T & force_magnitude ,
										   T & energy )
{
	effective_spring_constant = static_cast<T>(0) ;
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T force_component ;
	T delta_length_power_i = delta_length ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; ++i ) {
		force_component = force_length_parameters[i] * delta_length_power_i ;
		force_magnitude += force_component ;
		energy += force_component * delta_length / static_cast<T>(i+2) ;
		effective_spring_constant += force_component / ( static_cast<T>(i+1) * delta_length ) ;
		delta_length_power_i *= delta_length ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████  ██████  ██████   ██████ ███████         ███████ ██   ██ ██████   ██████  ███    ██ ███████ ███    ██ ████████ ██  █████  ██      
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██      ██    ██ ██   ██ ██      ██              ██       ██ ██  ██   ██ ██    ██ ████   ██ ██      ████   ██    ██    ██ ██   ██ ██      
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         █████   ██    ██ ██████  ██      █████           █████     ███   ██████  ██    ██ ██ ██  ██ █████   ██ ██  ██    ██    ██ ███████ ██      
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██      ██    ██ ██   ██ ██      ██              ██       ██ ██  ██      ██    ██ ██  ██ ██ ██      ██  ██ ██    ██    ██ ██   ██ ██      
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ██       ██████  ██   ██  ██████ ███████ ███████ ███████ ██   ██ ██       ██████  ██   ████ ███████ ██   ████    ██    ██ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::spring_force_exponential( const T & delta_length ,
										    const std::size_t & num_force_length_parameters ,
										    const std::vector<T> & force_length_parameters ,
											T & effective_spring_constant ,
										    T & force_magnitude ,
										    T & energy )
{
	effective_spring_constant = static_cast<T>(0) ;
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T force_component ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; i+=2 ) {
		force_component = force_length_parameters[i] * std::expm1( delta_length * force_length_parameters[i+1] ) ;
		force_magnitude += force_component ;
		energy += (force_component/force_length_parameters[i+1]) + (delta_length*force_length_parameters[i]) ;
		effective_spring_constant += force_length_parameters[i+1] * ( force_component + force_length_parameters[i] ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████  ██████  ██████   ██████ ███████         ██████   ██████  ██     ██ ███████ ██████  ██       █████  ██     ██ 
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██      ██    ██ ██   ██ ██      ██              ██   ██ ██    ██ ██     ██ ██      ██   ██ ██      ██   ██ ██     ██ 
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         █████   ██    ██ ██████  ██      █████           ██████  ██    ██ ██  █  ██ █████   ██████  ██      ███████ ██  █  ██ 
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██      ██    ██ ██   ██ ██      ██              ██      ██    ██ ██ ███ ██ ██      ██   ██ ██      ██   ██ ██ ███ ██ 
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ██       ██████  ██   ██  ██████ ███████ ███████ ██       ██████   ███ ███  ███████ ██   ██ ███████ ██   ██  ███ ███  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::spring_force_powerlaw( const T & delta_length ,
										 const std::size_t & num_force_length_parameters ,
										 const std::vector<T> & force_length_parameters ,
										 T & effective_spring_constant ,
										 T & force_magnitude ,
										 T & energy )
{
	effective_spring_constant = static_cast<T>(0) ;
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T force_component ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; i+=2 ) {
		force_component = force_length_parameters[i] * std::pow( delta_length , force_length_parameters[i+1] ) ;
		force_magnitude += force_component ;
		energy += force_component * delta_length / (force_length_parameters[i+1]+static_cast<T>(1)) ;
		effective_spring_constant += force_component / ( force_length_parameters[i+1] * delta_length ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████ 
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██    
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██    
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██              ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██    
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T Spring<T,N>::spring_stiffness_rest( void )
{
	std::size_t * num_force_length_parameters ;
	std::vector<T> * force_length_parameters ;
	Spring<T,N>::ForceLengthRelationship force_length_type ;
	T effective_spring_constant_rest = static_cast<T>(0) ;
	if( num_force_length_parameters_tension > 0 ) {
		num_force_length_parameters = & num_force_length_parameters_tension ;
		force_length_parameters     = & force_length_parameters_tension ;
		force_length_type           =   force_length_type_tension ;
	} else if( num_force_length_parameters_compression > 0 ) {
		num_force_length_parameters = & num_force_length_parameters_compression ;
		force_length_parameters     = & force_length_parameters_compression ;
		force_length_type           =   force_length_type_compression ;
	} else {
		return static_cast<T>(0) ;
	}
	switch( force_length_type ) {
		case Spring<T,N>::ForceLengthRelationship::polynomial :
			effective_spring_constant_rest = Spring<T,N>::spring_stiffness_rest_polynomial( * num_force_length_parameters ,
																							* force_length_parameters ) ;
			break ;
		case Spring<T,N>::ForceLengthRelationship::exponential :
			effective_spring_constant_rest = Spring<T,N>::spring_stiffness_rest_exponential( * num_force_length_parameters ,
																							 * force_length_parameters ) ;
			break ;
		case Spring<T,N>::ForceLengthRelationship::powerlaw :
			effective_spring_constant_rest = Spring<T,N>::spring_stiffness_rest_powerlaw( * num_force_length_parameters ,
																						  * force_length_parameters ) ;
			break ;
		default:
			effective_spring_constant_rest = static_cast<T>(0) ;
	}
	return effective_spring_constant_rest ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████         ██████   ██████  ██      ██    ██ ███    ██  ██████  ███    ███ ██  █████  ██      
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██            ██   ██ ██    ██ ██       ██  ██  ████   ██ ██    ██ ████  ████ ██ ██   ██ ██      
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██            ██████  ██    ██ ██        ████   ██ ██  ██ ██    ██ ██ ████ ██ ██ ███████ ██      
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██              ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██            ██      ██    ██ ██         ██    ██  ██ ██ ██    ██ ██  ██  ██ ██ ██   ██ ██      
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    ███████ ██       ██████  ███████    ██    ██   ████  ██████  ██      ██ ██ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T Spring<T,N>::spring_stiffness_rest_polynomial( const std::size_t & num_force_length_parameters ,
												 const std::vector<T> & force_length_parameters )
{
	return force_length_parameters[0] ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████         ███████ ██   ██ ██████   ██████  ███    ██ ███████ ███    ██ ████████ ██  █████  ██      
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██            ██       ██ ██  ██   ██ ██    ██ ████   ██ ██      ████   ██    ██    ██ ██   ██ ██      
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██            █████     ███   ██████  ██    ██ ██ ██  ██ █████   ██ ██  ██    ██    ██ ███████ ██      
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██              ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██            ██       ██ ██  ██      ██    ██ ██  ██ ██ ██      ██  ██ ██    ██    ██ ██   ██ ██      
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    ███████ ███████ ██   ██ ██       ██████  ██   ████ ███████ ██   ████    ██    ██ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T Spring<T,N>::spring_stiffness_rest_exponential( const std::size_t & num_force_length_parameters ,
												  const std::vector<T> & force_length_parameters )
{
	T spring_constant_rest = static_cast<T>(0) ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; i+=2 ) {
		spring_constant_rest += force_length_parameters[i] * force_length_parameters[i+1] ;
	}
	return spring_constant_rest ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████         ██████   ██████  ██     ██ ███████ ██████  ██       █████  ██     ██ 
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██            ██   ██ ██    ██ ██     ██ ██      ██   ██ ██      ██   ██ ██     ██ 
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██            ██████  ██    ██ ██  █  ██ █████   ██████  ██      ███████ ██  █  ██ 
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██              ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██            ██      ██    ██ ██ ███ ██ ██      ██   ██ ██      ██   ██ ██ ███ ██ 
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    ███████ ██       ██████   ███ ███  ███████ ██   ██ ███████ ██   ██  ███ ███  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
 T Spring<T,N>::spring_stiffness_rest_powerlaw( const std::size_t & num_force_length_parameters ,
											   const std::vector<T> & force_length_parameters )
{
	T spring_constant_rest = static_cast<T>(0) ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; i+=2 ) {
		if( force_length_parameters[i+1] == static_cast<T>(1.0) ) {
			spring_constant_rest += force_length_parameters[i] ;
		} else if( force_length_parameters[i+1] < static_cast<T>(1.0) ) {
			spring_constant_rest = static_cast<T>(1.0e+20) ;
			break ;
		}
	}
	return spring_constant_rest ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████           ██████  ███████ ████████         ███████  ██████  █████  ██      ███████         ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████ 
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██       ██         ██            ██      ██      ██   ██ ██      ██              ██         ██    ██ ██      ██      ████   ██ ██      ██      ██      
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ██   ███ █████      ██            ███████ ██      ███████ ██      █████           ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████ 
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██    ██ ██         ██                 ██ ██      ██   ██ ██      ██                   ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██ 
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████  ██████  ███████    ██    ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T Spring<T,N>::spring_get_scale_stiffness( void )
{
	T spring_scale_stiffness = effective_spring_constant ;
	if( spring_scale_stiffness < std::numeric_limits<T>::epsilon() ) {
		if( length >= rest_length ) {
			switch( force_length_type_tension ) {
				case Spring<T,N>::ForceLengthRelationship::polynomial :
					spring_scale_stiffness = force_length_parameters_tension[0] ;
					break ;
				case Spring<T,N>::ForceLengthRelationship::exponential :
					spring_scale_stiffness = force_length_parameters_tension[0] ;
					break ;
				case Spring<T,N>::ForceLengthRelationship::powerlaw :
					spring_scale_stiffness = force_length_parameters_tension[0] ;
					break ;
				default:
					break ;
			}
		} else {
			switch( force_length_type_compression ) {
				case Spring<T,N>::ForceLengthRelationship::polynomial :
					spring_scale_stiffness = force_length_parameters_compression[0] ;
					break ;
				case Spring<T,N>::ForceLengthRelationship::exponential :
					spring_scale_stiffness = force_length_parameters_compression[0] ;
					break ;
				case Spring<T,N>::ForceLengthRelationship::powerlaw :
					spring_scale_stiffness = force_length_parameters_compression[0] ;
					break ;
				default:
					break ;
			}
		}
	}
	return spring_scale_stiffness ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ██████  ██████  ██ ███    ██  ██████          ██████  ███████ ███████  ██████  █████  ██      ███████         ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████ 
//  ██      ██   ██ ██   ██ ██ ████   ██ ██               ██   ██ ██      ██      ██      ██   ██ ██      ██              ██         ██    ██ ██      ██      ████   ██ ██      ██      ██      
//  ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ██████  █████   ███████ ██      ███████ ██      █████           ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████ 
//       ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██   ██ ██           ██ ██      ██   ██ ██      ██                   ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██ 
//  ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ██   ██ ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::spring_rescale( const T & scale_length , const T & scale_stiffness )
{
	rest_length *= scale_length ;
	switch( force_length_type_tension ) {
		case Spring<T,N>::ForceLengthRelationship::polynomial :
			for( std::size_t i = 0 ; i < num_force_length_parameters_tension ; ++i ) {
				force_length_parameters_tension[i] *= scale_stiffness*std::pow(scale_length,static_cast<T>(1)-static_cast<T>(i)) ; ;
			}
			break ;
		case Spring<T,N>::ForceLengthRelationship::exponential :
			for( std::size_t i = 0 ; i < num_force_length_parameters_tension ; i+=2 ) {
				force_length_parameters_tension[i] *= scale_stiffness*scale_length ;
				force_length_parameters_tension[i+1] /= scale_length ;
			}
			break ;
		case Spring<T,N>::ForceLengthRelationship::powerlaw :
			for( std::size_t i = 0 ; i < num_force_length_parameters_tension ; i+=2 ) {
				force_length_parameters_tension[i] *= scale_stiffness*std::pow(scale_length,static_cast<T>(1)-force_length_parameters_tension[i+1]) ;
			}
			break ;
		default:
			break ;
	}
	switch( force_length_type_compression ) {
		case Spring<T,N>::ForceLengthRelationship::polynomial :
			for( std::size_t i = 0 ; i < num_force_length_parameters_compression ; ++i ) {
				force_length_parameters_compression[i] *= scale_stiffness*std::pow(scale_length,static_cast<T>(1)-static_cast<T>(i)) ; ;
			}
			break ;
		case Spring<T,N>::ForceLengthRelationship::exponential :
			for( std::size_t i = 0 ; i < num_force_length_parameters_compression ; i+=2 ) {
				force_length_parameters_compression[i] *= scale_stiffness*scale_length ;
				force_length_parameters_compression[i+1] /= scale_length ;
			}
			break ;
		case Spring<T,N>::ForceLengthRelationship::powerlaw :
			for( std::size_t i = 0 ; i < num_force_length_parameters_compression ; i+=2 ) {
				force_length_parameters_compression[i] *= scale_stiffness*std::pow(scale_length,static_cast<T>(1)-force_length_parameters_compression[i+1]) ;
			}
			break ;
		default:
			break ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   ad88888ba                            88                            888b      88                                                                    88         
//  d8"     "8b                           ""                            8888b     88                ,d                                                  88         
//  Y8,                                                                 88 `8b    88                88                                                  88         
//  `Y8aaaaa,    8b,dPPYba,   8b,dPPYba,  88  8b,dPPYba,    ,adPPYb,d8  88  `8b   88   ,adPPYba,  MM88MMM  8b      db      d8   ,adPPYba,   8b,dPPYba,  88   ,d8   
//    `"""""8b,  88P'    "8a  88P'   "Y8  88  88P'   `"8a  a8"    `Y88  88   `8b  88  a8P_____88    88     `8b    d88b    d8'  a8"     "8a  88P'   "Y8  88 ,a8"    
//          `8b  88       d8  88          88  88       88  8b       88  88    `8b 88  8PP"""""""    88      `8b  d8'`8b  d8'   8b       d8  88          8888[      
//  Y8a     a8P  88b,   ,a8"  88          88  88       88  "8a,   ,d88  88     `8888  "8b,   ,aa    88,      `8bd8'  `8bd8'    "8a,   ,a8"  88          88`"Yba,   
//   "Y88888P"   88`YbbdP"'   88          88  88       88   `"YbbdP"Y8  88      `888   `"Ybbd8"'    "Y888      YP      YP       `"YbbdP"'   88          88   `Y8a  
//               88                                         aa,    ,88                                                                                             
//               88                                          "Y8bbdP"                                                                                              
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   ██████ ██████  ███████  █████  ████████ ███████         ███████ ██████  ██████  ██ ███    ██  ██████          ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██          ██████  ██████       ██ 
//  ██      ██   ██ ██      ██   ██    ██    ██              ██      ██   ██ ██   ██ ██ ████   ██ ██               ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██          ██    ██ ██   ██      ██ 
//  ██      ██████  █████   ███████    ██    █████           ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████           ██    ██ ██████       ██ 
//  ██      ██   ██ ██      ██   ██    ██    ██                   ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██          ██    ██ ██   ██ ██   ██ 
//   ██████ ██   ██ ███████ ██   ██    ██    ███████ ███████ ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ ███████  ██████  ██████   █████  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<ASpringNetwork> ASpringNetwork::create_spring_network_obj( const NetworkParameters & network_parameters )
{
	std::unique_ptr<ASpringNetwork> asn = NULL ;
	std::cout
		<< '\t' << network_parameters.precision
		<< '\t' << network_parameters.num_dimensions << 'D'
		<< std::endl ;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ███████ ████████ ██    ██ ██████  
//  ██      ██         ██    ██    ██ ██   ██ 
//  ███████ █████      ██    ██    ██ ██████  
//       ██ ██         ██    ██    ██ ██      
//  ███████ ███████    ██     ██████  ██      
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
	uni_0_1 = std::uniform_real_distribution<T>( 0.0,+1.0) ;
	uni_1_1 = std::uniform_real_distribution<T>(-1.0,+1.0) ;
	num_points  = network_parameters.num_points ;
	num_springs = network_parameters.num_springs ;
	points  = std::vector<  Point<T,N> >( num_points  ) ;
	springs = std::vector< Spring<T,N> >( num_springs ) ;
	nodes   = std::vector<        Node >( num_points  ) ;
	load_network_binary( file_input_nodes.c_str() , file_input_springs.c_str() ) ;
	construct_network() ;
	//
	if( network_parameters.parallelism_enabled ) {
		if( network_parameters.user_num_threads > 1 ) {
			if( num_springs > 1e4 ) {
				num_threads = std::min( static_cast<int>(network_parameters.user_num_threads) , omp_get_max_threads() ) ;
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
		#ifdef _WIN32
			putenv( "OMP_WAIT_POLICY=ACTIVE" ) ;
			putenv( "OMP_DYNAMIC=FALSE" ) ;
			putenv( "OMP_PROC_BIND=TRUE" ) ;
		#else
			setenv( "OMP_WAIT_POLICY" , "ACTIVE" , 1 ) ;
			setenv( "OMP_DYNAMIC" , "FALSE" , 1 ) ;
			setenv( "OMP_PROC_BIND" , "TRUE" , 1 ) ;
		#endif
	}
	//
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██       ██████   █████  ██████          ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██         ██████  ██ ███    ██  █████  ██████  ██    ██ 
//  ██      ██    ██ ██   ██ ██   ██         ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ████   ██ ██   ██ ██   ██  ██  ██  
//  ██      ██    ██ ███████ ██   ██         ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████           ██████  ██ ██ ██  ██ ███████ ██████    ████   
//  ██      ██    ██ ██   ██ ██   ██         ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ██  ██ ██ ██   ██ ██   ██    ██    
//  ███████  ██████  ██   ██ ██████  ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ ███████ ██████  ██ ██   ████ ██   ██ ██   ██    ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::load_network_binary( const char * file_nodes ,
											  const char * file_springs )
{
	std::FILE * file_ptr ;
	T             * read_data_T ;
	std::uint32_t * read_data_uint32 ;
	std::uint8_t  * read_data_uint8 ;
	std::size_t NFLPT ;
	std::size_t NFLPC ;
	
	// points
	if( file_nodes != NULL ) {
		file_ptr = std::fopen( file_nodes , "rb" ) ;
		read_data_T     = new T            [N+N] ;
		read_data_uint8 = new std::uint8_t [N  ] ;
		for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
			std::fread( read_data_T     , sizeof(T)            , N+N , file_ptr ) ;
			std::fread( read_data_uint8 , sizeof(std::uint8_t) , N   , file_ptr ) ;
			p->fixed_all_dim = true ;
			for( std::size_t n = 0 ; n < N ; ++n ) {
				p->position[n]      = read_data_T[n  ] ;
				p->force_applied[n] = read_data_T[n+N] ;
				p->fixed_dim[n] = ( read_data_uint8[n] != 0 ) ;
				p->fixed_all_dim = p->fixed_all_dim && p->fixed_dim[n] ;
			}
		}
		delete [] read_data_T     ;
		delete [] read_data_uint8 ;
		std::fclose( file_ptr ) ;
	}
	
	// springs
	if( file_springs != NULL ) {
		file_ptr = std::fopen( file_springs , "rb" ) ;
		read_data_uint32  = new std::uint32_t [4] ;
		read_data_uint8   = new std::uint8_t  [2] ;
		for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
			//
			std::fread( &read_data_uint32[0] , sizeof(std::uint32_t) , 2 , file_ptr ) ;
			s->start = &points[ read_data_uint32[0] ] ;
			s->end   = &points[ read_data_uint32[1] ] ;
			//
			std::fread(  read_data_uint8     , sizeof(std::uint8_t)  , 2 , file_ptr ) ;
			s->force_length_type_tension     = static_cast<typename Spring<T,N>::ForceLengthRelationship>(read_data_uint8[0]) ;
			s->force_length_type_compression = static_cast<typename Spring<T,N>::ForceLengthRelationship>(read_data_uint8[1]) ;
			//
			read_data_T = new T [1] ;
			std::fread( read_data_T      , sizeof(T) , 1 , file_ptr ) ;
			s->rest_length = read_data_T[0] ;
			delete [] read_data_T ;
			//
			std::fread( &read_data_uint32[2] , sizeof(std::uint32_t) , 1 , file_ptr ) ;
			NFLPT = static_cast<std::size_t>(read_data_uint32[2]) ;
			read_data_T = new T [NFLPT] ;
			std::fread( read_data_T , sizeof(T) , NFLPT , file_ptr ) ;
			s->num_force_length_parameters_tension     = NFLPT ;
			s->force_length_parameters_tension = std::vector<T>( &read_data_T[0] , &read_data_T[NFLPT] ) ;
			delete [] read_data_T ;
			//
			std::fread( &read_data_uint32[3]  , sizeof(std::uint32_t) , 1 , file_ptr ) ;
			NFLPC = static_cast<std::size_t>(read_data_uint32[3]) ;
			read_data_T = new T [NFLPC] ;
			std::fread( read_data_T , sizeof(T) , NFLPC , file_ptr ) ;
			s->num_force_length_parameters_compression = NFLPC ;
			s->force_length_parameters_compression = std::vector<T>( &read_data_T[0] , &read_data_T[NFLPC] ) ;
			delete [] read_data_T ;
		}
		delete [] read_data_uint32 ;
		delete [] read_data_uint8  ;
		std::fclose( file_ptr ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████  █████  ██    ██ ███████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██         ██████  ██ ███    ██  █████  ██████  ██    ██ 
//  ██      ██   ██ ██    ██ ██              ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ████   ██ ██   ██ ██   ██  ██  ██  
//  ███████ ███████ ██    ██ █████           ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████           ██████  ██ ██ ██  ██ ███████ ██████    ████   
//       ██ ██   ██  ██  ██  ██              ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ██  ██ ██ ██   ██ ██   ██    ██    
//  ███████ ██   ██   ████   ███████ ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ ███████ ██████  ██ ██   ████ ██   ██ ██   ██    ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::save_network_binary( const char * file_nodes ,
											  const char * file_springs )
{
	std::FILE * file_ptr ;
	T             * write_data_T ;
	std::uint32_t * write_data_uint32 ;
	std::uint8_t  * write_data_uint8 ;
	std::size_t NFLPT ;
	std::size_t NFLPC ;
	
	// points
	if( file_nodes != NULL ) {
		file_ptr = std::fopen( file_nodes , "wb" ) ;
		write_data_T     = new T            [N+N] ;
		write_data_uint8 = new std::uint8_t [N  ] ;
		for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
			for( std::size_t n = 0 ; n < N ; ++n ) {
				write_data_T[n  ] = p->position[n]      ;
				write_data_T[n+N] = p->force_applied[n] ;
				write_data_uint8[n] = static_cast<std::uint8_t>( p->fixed_dim[n] ) ;
			}
			std::fwrite( write_data_T     , sizeof(T)            , N+N , file_ptr ) ;
			std::fwrite( write_data_uint8 , sizeof(std::uint8_t) , N   , file_ptr ) ;
		}
		delete [] write_data_T     ;
		delete [] write_data_uint8 ;
		std::fclose( file_ptr ) ;
	}
	
	// springs
	if( file_springs != NULL ) {
		file_ptr = std::fopen( file_springs , "wb" ) ;
		write_data_uint32 = new std::uint32_t [4] ;
		write_data_uint8  = new std::uint8_t  [2] ;
		for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
			NFLPT = s->num_force_length_parameters_tension ;
			NFLPC = s->num_force_length_parameters_compression ;
			write_data_T = new T [1+NFLPT+NFLPC] ;
			write_data_uint32[0] = static_cast<std::uint32_t>( s->start - &points[0] ) ;
			write_data_uint32[1] = static_cast<std::uint32_t>( s->end   - &points[0] ) ;
			write_data_uint32[2] = static_cast<std::uint32_t>( NFLPT ) ;
			write_data_uint32[3] = static_cast<std::uint32_t>( NFLPC ) ;
			write_data_uint8[0] = static_cast<std::uint8_t>( s->force_length_type_tension     ) ;
			write_data_uint8[1] = static_cast<std::uint8_t>( s->force_length_type_compression ) ;
			write_data_T[0] = s->rest_length ;
			std::copy( s->force_length_parameters_tension.begin()     , s->force_length_parameters_tension.end()     , &write_data_T[1      ] ) ;
			std::copy( s->force_length_parameters_compression.begin() , s->force_length_parameters_compression.end() , &write_data_T[1+NFLPT] ) ;
			std::fwrite( &write_data_uint32[0]  , sizeof(std::uint32_t) , 2     , file_ptr ) ;
			std::fwrite(  write_data_uint8      , sizeof(std::uint8_t)  , 2     , file_ptr ) ;
			std::fwrite( &write_data_T[0]       , sizeof(T)             , 1     , file_ptr ) ;
			std::fwrite( &write_data_uint32[2]  , sizeof(std::uint32_t) , 1     , file_ptr ) ;
			std::fwrite( &write_data_T[1]       , sizeof(T)             , NFLPT , file_ptr ) ;
			std::fwrite( &write_data_uint32[3]  , sizeof(std::uint32_t) , 1     , file_ptr ) ;
			std::fwrite( &write_data_T[1+NFLPT] , sizeof(T)             , NFLPC , file_ptr ) ;
			delete [] write_data_T ;
		}
		delete [] write_data_uint32 ;
		delete [] write_data_uint8  ;
		std::fclose( file_ptr ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██ 
// ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██            ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██  
// ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██            ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████   
// ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██            ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██  
//  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██    ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::construct_network( void )
{
	for( std::size_t p = 0 ; p < num_points ; ++p ) {
		nodes[p].point = &points[p] ;
		points[p].not_referenced = true ;
	}
	std::size_t p_start ;
	std::size_t p_end ;
	for( std::size_t s = 0 ; s < num_springs ; ++s ) {
		p_start = springs[s].start - &points[0] ;
		p_end   = springs[s].end   - &points[0] ;
		// ignore link if spring has no force-length relationships, or if both points are fixed in all dimensions
		bool both_points_fixed = points[p_start].fixed_all_dim
		                      && points[p_end  ].fixed_all_dim ;
		bool none_force_length = springs[s].force_length_type_tension    ==Spring<T,N>::ForceLengthRelationship::none
		                      && springs[s].force_length_type_compression==Spring<T,N>::ForceLengthRelationship::none ;
		if( !both_points_fixed && !none_force_length ) {
			nodes[p_start].links.push_back( (Link){ -1 , &points[p_end  ] , &springs[s] , static_cast<T>(+1) } ) ;
			nodes[p_start].links.push_back( (Link){ -1 , &points[p_start] , &springs[s] , static_cast<T>(-1) } ) ;
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
// ███████  ██████  ██      ██    ██ ███████ 
// ██      ██    ██ ██      ██    ██ ██      
// ███████ ██    ██ ██      ██    ██ █████   
//      ██ ██    ██ ██       ██  ██  ██      
// ███████  ██████  ███████   ████   ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::solve( void )
{
	get_scale_network() ;
	rescale_network(+1) ;
	if( algorithm=="newton"   ) { minimize_energy_newton() ; } else
	if( algorithm=="steepest" ) { minimize_energy()        ; } else
	if( algorithm=="anneal"   ) { anneal()                 ; }
	rescale_network(-1) ;
	save_output() ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ███████ ████████         ███████  ██████  █████  ██      ███████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██ 
// ██       ██         ██            ██      ██      ██   ██ ██      ██              ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██  
// ██   ███ █████      ██            ███████ ██      ███████ ██      █████           ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████   
// ██    ██ ██         ██                 ██ ██      ██   ██ ██      ██              ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██  
//  ██████  ███████    ██    ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::get_scale_network( void )
{
	// want to reset minimum positions to the origin
	scale_position_min = nodes.begin()->point->position ;
	scale_position_max = nodes.begin()->point->position ;
	iterNode n ;
	iterNode n_end ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( std::size_t d = 0 ; d < N ; ++d ) {
			if( n->point->position[d] < scale_position_min[d] ) {
				scale_position_min[d] = n->point->position[d] ;
			} else if( n->point->position[d] > scale_position_max[d] ) {
				scale_position_max[d] = n->point->position[d] ;
			}
		}
	}
	scale_position_range = scale_position_max - scale_position_min ;

	// want to scale network size and stiffnesses by average spring length and stiffness
	scale_stiffness = static_cast<T>(0.0) ;
	scale_length = static_cast<T>(0.0) ;
	iterSpring s ;
	iterSpring s_end ;
	for( s = springs.begin() , s_end = springs.end() ; s != s_end ; ++s ) {
		s->spring_energy() ;
		scale_stiffness += s->spring_get_scale_stiffness() ;
		scale_length += s->length ;
	}
	scale_stiffness /= static_cast<T>(springs.size()) ;
	scale_length /= static_cast<T>(springs.size()) ;
	scale_force = scale_length * scale_stiffness ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██████  ███████ ███████  ██████  █████  ██      ███████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██ 
// ██   ██ ██      ██      ██      ██   ██ ██      ██              ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██  
// ██████  █████   ███████ ██      ███████ ██      █████           ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████   
// ██   ██ ██           ██ ██      ██   ██ ██      ██              ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██  
// ██   ██ ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::rescale_network( const int & direction )
{
	T scale_length_inv ;
	T scale_force_inv ;
	T scale_stiffness_inv ;
	switch( direction ) {
		case +1:
			// set to normalize scale from original scale
			scale_length_inv = static_cast<T>(1) / scale_length ;
			scale_force_inv = static_cast<T>(1) / scale_force ;
			scale_stiffness_inv = static_cast<T>(1) / scale_stiffness ;
			for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
				p->position -= scale_position_min ;
				p->position *= scale_length_inv ;
				p->force_applied *= scale_force_inv ;
			}
			for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
				s->spring_rescale( scale_length_inv , scale_stiffness_inv ) ;
			}
			break ;

		case -1:
			// reset to original scale from normalize scale
			for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
				p->position *= scale_length ;
				p->position += scale_position_min ;
				p->force_applied *= scale_force ;
			}
			for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
				s->spring_rescale( scale_length , scale_stiffness ) ;
			}
			break ;

		default:
			break ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███████ ██ ███    ██ ██████          ███    ███  █████  ██   ██         ███████ ██████  ██████  ██ ███    ██  ██████          ██      ███████ ███    ██  ██████  ████████ ██   ██ 
// ██      ██ ████   ██ ██   ██         ████  ████ ██   ██  ██ ██          ██      ██   ██ ██   ██ ██ ████   ██ ██               ██      ██      ████   ██ ██          ██    ██   ██ 
// █████   ██ ██ ██  ██ ██   ██         ██ ████ ██ ███████   ███           ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ██      █████   ██ ██  ██ ██   ███    ██    ███████ 
// ██      ██ ██  ██ ██ ██   ██         ██  ██  ██ ██   ██  ██ ██               ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██      ██      ██  ██ ██ ██    ██    ██    ██   ██ 
// ██      ██ ██   ████ ██████  ███████ ██      ██ ██   ██ ██   ██ ███████ ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ███████ ███████ ██   ████  ██████     ██    ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::find_max_spring_length( void )
{
	Vector<T,N> x_min = points.begin()->position ;
	Vector<T,N> x_max = points.begin()->position ;
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
	max_spring_length *= 3.0 ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ████████  ██████  ████████  █████  ██              ███████ ███    ██ ███████ ██████   ██████  ██    ██ 
//    ██    ██    ██    ██    ██   ██ ██              ██      ████   ██ ██      ██   ██ ██        ██  ██  
//    ██    ██    ██    ██    ███████ ██              █████   ██ ██  ██ █████   ██████  ██   ███   ████   
//    ██    ██    ██    ██    ██   ██ ██              ██      ██  ██ ██ ██      ██   ██ ██    ██    ██    
//    ██     ██████     ██    ██   ██ ███████ ███████ ███████ ██   ████ ███████ ██   ██  ██████     ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T SpringNetwork<T,N>::total_energy( void )
{
	// using klein summation, modification of kahan summation,
	// to account for floating point precision issues with sums
	// of many small numbers, or small + large numbers
	ksum.reset() ;
	iterNode n ;
	iterNode n_end ;
	iterLink l ;
	iterLink l_end ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			ksum.add( l->spring->energy ) ;
		}
	}
	/*
	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
		ksum.add( s->energy ) ;
	}
	//*/
	T energy = ksum.result() ;
	/*
	// external work done by applied loads
	Vector<T,N> displacement ;
	for( std::size_t p = 0 ; p < num_points ; ++p ) {
		displacement = points[p].position ;
		displacement -= points_init[p].position ;
		energy -= points[p].force_applied.dot( displacement ) ;
	}
	//*/
	return energy ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██    ██ ██████  ██████   █████  ████████ ███████         ███████ ██████  ██████  ██ ███    ██  ██████  ███████ 
// ██    ██ ██   ██ ██   ██ ██   ██    ██    ██              ██      ██   ██ ██   ██ ██ ████   ██ ██       ██      
// ██    ██ ██████  ██   ██ ███████    ██    █████           ███████ ██████  ██████  ██ ██ ██  ██ ██   ███ ███████ 
// ██    ██ ██      ██   ██ ██   ██    ██    ██                   ██ ██      ██   ██ ██ ██  ██ ██ ██    ██      ██ 
//  ██████  ██      ██████  ██   ██    ██    ███████ ███████ ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::update_springs( void )
{
	// internal energy due to spring tension/compression
	// parallelize computation for more than 10^4 springs
	#pragma omp for schedule(static)
	for( std::size_t i = 0 ; i < springs_used.size() ; ++i ) {
		springs_used[i]->spring_energy() ;
	}
	/*
	#pragma omp for schedule(static)
	for( std::size_t s = 0 ; s < num_springs ; ++s ) {
		springs[s].spring_energy() ;
	}
	//*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██    ██ ██████  ██████   █████  ████████ ███████         ███████  ██████  ██████   ██████ ███████ ███████ 
// ██    ██ ██   ██ ██   ██ ██   ██    ██    ██              ██      ██    ██ ██   ██ ██      ██      ██      
// ██    ██ ██████  ██   ██ ███████    ██    █████           █████   ██    ██ ██████  ██      █████   ███████ 
// ██    ██ ██      ██   ██ ██   ██    ██    ██              ██      ██    ██ ██   ██ ██      ██           ██ 
//  ██████  ██      ██████  ██   ██    ██    ███████ ███████ ██       ██████  ██   ██  ██████ ███████ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::update_forces( void )
{
	/*
	// account for force imbalances, "potential energy"
	// combine external applied force & internal spring forces
	iterNode n ;
	iterNode n_end ;
	iterLink l ;
	iterLink l_end ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		n->point->net_force = n->point->force_applied ;
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			if( l->spring_direction > static_cast<T>(0) ) {
				n->point->net_force += l->spring->get_force() ;
			} else {
				n->point->net_force -= l->spring->get_force() ;
			}
		}
		n->point->net_force_magnitude = n->point->net_force.norm() ;
	}
	sum_net_force_magnitude = static_cast<T>(0) ;
	for( iterNode n = nodes.begin() ; n != nodes.end() ; ++n ) {
		sum_net_force_magnitude += n->point->net_force_magnitude ;
	}
	//*/

	// net force on each point (including fixed points!)
	/*
	#pragma omp for schedule(static)
	for( std::size_t p = 0 ; p < num_points ; ++p ) {
		points[p].net_force = static_cast<T>(0.0) ;
	}
	//*/
	#pragma omp single
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		p->net_force = static_cast<T>(0.0) ;
	}
	#pragma omp single
	for( std::size_t i = 0 ; i < springs_used.size() ; ++i ) {
		springs_used[i]->start->net_force += springs_used[i]->force ;
		springs_used[i]->end  ->net_force -= springs_used[i]->force ;
	}
	/*
	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
		s->start->net_force += s->force ;
		s->end->net_force   -= s->force ;
	}
	//*/
	#pragma omp for schedule(static)
	for( std::size_t p = 0 ; p < num_points ; ++p ) {
		Point<T,N> * current_point = &points[p] ;
		current_point->net_force += current_point->force_applied ;
		current_point->net_force_magnitude = current_point->net_force.norm() ;
	}
	/*
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		p->net_force = static_cast<T>(0.0) ;
	}
	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
		Vector<T,N> spring_force = s->get_force() ;
		s->start->net_force += spring_force ;
		s->end->net_force   -= spring_force ;
	}
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		p->net_force += p->force_applied ;
		p->net_force_magnitude = p->net_force.norm() ;
	}
	//*/
	#pragma omp single
	{
		ksum.reset() ;
		max_net_force_magnitude = static_cast<T>(0.0) ;
		if( include_force_fixed_nodes ) {
			for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
				if( !p->not_referenced ) {
					max_net_force_magnitude = ( p->net_force_magnitude > max_net_force_magnitude ) ? p->net_force_magnitude : max_net_force_magnitude ;
					ksum.add( p->net_force_magnitude ) ;
				}
			}
		} else {
			for( iterNode n = nodes.begin() ; n != nodes.end() ; ++n ) {
				max_net_force_magnitude = ( n->point->net_force_magnitude > max_net_force_magnitude ) ? n->point->net_force_magnitude : max_net_force_magnitude ;
				ksum.add( n->point->net_force_magnitude ) ;
			}
		}
		sum_net_force_magnitude = ksum.result() ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ███████ ████████          ██████  ██████       ██ ███████  ██████ ████████ ██ ██    ██ ███████ 
// ██       ██         ██            ██    ██ ██   ██      ██ ██      ██         ██    ██ ██    ██ ██      
// ██   ███ █████      ██            ██    ██ ██████       ██ █████   ██         ██    ██ ██    ██ █████   
// ██    ██ ██         ██            ██    ██ ██   ██ ██   ██ ██      ██         ██    ██  ██  ██  ██      
//  ██████  ███████    ██    ███████  ██████  ██████   █████  ███████  ██████    ██    ██   ████   ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T SpringNetwork<T,N>::get_objective( void )
{
	if( objective=="energy"   ) { return total_energy()          ; } else
	if( objective=="sumforce" ) { return sum_net_force_magnitude ; } else
	if( objective=="maxforce" ) { return max_net_force_magnitude ; }
	else { return static_cast<T>(0.0) ; }
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
// ███    ███  ██████  ██    ██ ███████         ██████   ██████  ██ ███    ██ ████████ ███████         ███████  ██████  ██████   ██████ ███████ 
// ████  ████ ██    ██ ██    ██ ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ██      ██    ██ ██   ██ ██      ██      
// ██ ████ ██ ██    ██ ██    ██ █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         █████   ██    ██ ██████  ██      █████   
// ██  ██  ██ ██    ██  ██  ██  ██              ██      ██    ██ ██ ██  ██ ██    ██         ██         ██      ██    ██ ██   ██ ██      ██      
// ██      ██  ██████    ████   ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ██       ██████  ██   ██  ██████ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_force( const T & step_size )
{
	// apply net force & move small displacement towards equilibrating position
	/*
	iterNode n ;
	iterNode n_end ;
	Vector<T,N> small_displacement ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		small_displacement = n->point->net_force ;
		small_displacement *= step_size ;
		n->point->move( small_displacement ) ;
	}
	//*/
	Vector<T,N> small_displacement ;
	#pragma omp for schedule(static) private(small_displacement)
	for( std::size_t n = 0 ; n < nodes.size() ; ++n ) {
		small_displacement = nodes[n].point->net_force ;
		small_displacement *= step_size ;
		nodes[n].point->move( small_displacement ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███    ███  ██████  ██    ██ ███████         ██████   ██████  ██ ███    ██ ████████ ███████         ██████   █████  ███    ██ ██████  
// ████  ████ ██    ██ ██    ██ ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ██   ██ ██   ██ ████   ██ ██   ██ 
// ██ ████ ██ ██    ██ ██    ██ █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         ██████  ███████ ██ ██  ██ ██   ██ 
// ██  ██  ██ ██    ██  ██  ██  ██              ██      ██    ██ ██ ██  ██ ██    ██         ██         ██   ██ ██   ██ ██  ██ ██ ██   ██ 
// ██      ██  ██████    ████   ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ██   ██ ██   ██ ██   ████ ██████  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_rand( const T & amplitude )
{
	iterNode n ;
	iterNode n_end ;
	Vector<T,N> small_displacement ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( std::size_t d = 0 ; d < N ; ++d ) {
			small_displacement[d] = amplitude * uni_1_1(rng) ;
		}
		n->point->move( small_displacement ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ████████ ███████ ███████ ████████         ██████  ███████ ██████   ██████   ██████  ████████ 
//    ██    ██      ██         ██            ██   ██ ██      ██   ██ ██    ██ ██    ██    ██    
//    ██    █████   ███████    ██            ██████  █████   ██████  ██    ██ ██    ██    ██    
//    ██    ██           ██    ██            ██   ██ ██      ██   ██ ██    ██ ██    ██    ██    
//    ██    ███████ ███████    ██    ███████ ██   ██ ███████ ██████   ██████   ██████     ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  █████   ██████  ██████ ███████ ██████  ████████         ███    ██ ███████ ██     ██         ██████   ██████  ██ ███    ██ ████████ ███████ 
// ██   ██ ██      ██      ██      ██   ██    ██            ████   ██ ██      ██     ██         ██   ██ ██    ██ ██ ████   ██    ██    ██      
// ███████ ██      ██      █████   ██████     ██            ██ ██  ██ █████   ██  █  ██         ██████  ██    ██ ██ ██ ██  ██    ██    ███████ 
// ██   ██ ██      ██      ██      ██         ██            ██  ██ ██ ██      ██ ███ ██         ██      ██    ██ ██ ██  ██ ██    ██         ██ 
// ██   ██  ██████  ██████ ███████ ██         ██    ███████ ██   ████ ███████  ███ ███  ███████ ██       ██████  ██ ██   ████    ██    ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██   ██ ███████  █████  ████████         ██    ██ ██████  
// ██   ██ ██      ██   ██    ██            ██    ██ ██   ██ 
// ███████ █████   ███████    ██            ██    ██ ██████  
// ██   ██ ██      ██   ██    ██            ██    ██ ██      
// ██   ██ ███████ ██   ██    ██    ███████  ██████  ██      
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T SpringNetwork<T,N>::heat_up( const T & amplitude , const std::size_t & num_config_test )
{
	T energy = static_cast<T>(0) ;
	T energy_max = std::numeric_limits<T>::min() ;
	points_prev = points ;
	for( std::size_t i = 0 ; i < num_config_test ; ++i ) {
		move_points_rand( amplitude ) ;
		update_springs() ;
		update_forces() ; // not needed if energy does not depend on forces
		energy = get_objective() ;
		if( energy > energy_max ) {
			energy_max = energy ;
		}
		points = points_prev ;
	}
	return energy_max ;
}

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
	#pragma omp parallel
	for( std::size_t iter = 0 ; iter < local_num_iter_max ; ++iter ) {

		#pragma omp single
		{
			// perturb the system with some randomness
			if( is_shaking ) {
				if( iter%num_iter_shake == 0 ) {
					move_points_rand( shake_step_size ) ;
					shake_step_size *= shake_step_size_reduction ;
					update_springs() ;
					update_forces() ; // also computes sum_net_force_magnitude
				}
			}
		}

		#pragma omp single
		{
			// basic line search
			energy_prev = energy ;
			points_prev = points ;
			step_size /= step_size_reduction ;
			step_size /= step_size_reduction ;
		}
		do {
			#pragma omp single
			{
				step_size *= step_size_reduction ;
				points = points_prev ;
			}
			move_points_force( step_size ) ;
			update_springs() ;
			update_forces() ; // also computes sum_net_force_magnitude
			#pragma omp single
			{
				energy = get_objective() ;
			}
		} while( (energy>energy_prev) && (step_size>step_size_min) ) ;

		#pragma omp single
		{
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
		}
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
	#pragma omp parallel
	for( std::size_t iter = 0 ; iter < local_num_iter_max ; ++iter ) {

		#pragma omp single
		{
			// basic line search
			energy_prev = energy ;
			points_prev = points ;
			step_size /= step_size_reduction ;
			step_size /= step_size_reduction ;
			compute_newton_step_direction() ;
		}
		do {
			#pragma omp single
			{
				step_size *= step_size_reduction ;
				step_size = ( step_size > step_size_max ) ? step_size_max : step_size ;
				points = points_prev ;
			}
			move_points_force( step_size ) ; // move_points_newton( step_size ) ;
			update_springs() ;
			update_forces() ; // also computes sum_net_force_magnitude
			#pragma omp single
			{
				energy = get_objective() ;
			}
		} while( (energy>=energy_prev) && (step_size>step_size_min) ) ;

		#pragma omp single
		{
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
		}
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
//  █████  ███    ██ ███    ██ ███████  █████  ██      
// ██   ██ ████   ██ ████   ██ ██      ██   ██ ██      
// ███████ ██ ██  ██ ██ ██  ██ █████   ███████ ██      
// ██   ██ ██  ██ ██ ██  ██ ██ ██      ██   ██ ██      
// ██   ██ ██   ████ ██   ████ ███████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::anneal( void )
{
	// copy current state to initial, previous, and best states
	points_init = points ;
	points_prev = points ;
	points_best = points ;
	update_springs() ;
	update_forces() ; // also computes sum_net_force_magnitude
	T energy = get_objective() ;
	//T energy_init = energy ;
	T energy_prev = energy ;
	T energy_best = energy ;

	//
	double time_start ;
	double time_prev ;
	double time_curr ;
	
	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 500000 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :      0 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :   5000 ;
	T local_tolerance_sum_net_force    = ( tolerance_sum_net_force    > 0.0 ) ? tolerance_sum_net_force    : 1E-12 ;
	T local_tolerance_change_objective = ( tolerance_change_objective > 0.0 ) ? tolerance_change_objective : 1E-12 ;

	// annealing solver parameters
	// TODO // get these parameters from NetworkParameters? (default/file)
	T step_size = static_cast<T>(0.2) ;
	T step_size_reduction = 0.9 ;
	T temperature_reduction = static_cast<T>(0.99) ;
	T energy_compare ;
	T energy_compare_min = static_cast<T>(1E-12) ;
	T relative_change_energy ;
	T mean_energy_change = static_cast<T>(0.0) ;
	T sum_net_force_magnitude_prev = sum_net_force_magnitude ;
	//std::size_t num_iter_heatup = local_num_iter_max / 5 ;
	std::size_t num_consecutive_reject     =   0 ;
	std::size_t num_consecutive_reject_max = 200 ;
	std::size_t num_since_last_best     = 0 ;
	std::size_t num_since_last_best_max = local_num_iter_max / 100 ;
	std::size_t num_small_change     =   0 ;
	std::size_t num_small_change_max = 200 ;
	std::size_t num_temperature_reductions     =   0 ;
	std::size_t num_temperature_reductions_max = 300 ;
	bool reboot ;
	bool break_temperature_reductions = false ;
	bool break_sum_net_force_magnitude = false ;
	bool break_consecutive_reject = false ;
	bool break_isnan_energy = false ;
	find_max_spring_length() ;
	
	// choose starting temperature
	// equivalent to 1/e probability to accept energy difference
	T heatup_amplitude = 0.01 ;
	std::size_t heatup_num_config = 1000 ;
	time_start = omp_get_wtime() ;
	energy_compare = heat_up( heatup_amplitude , heatup_num_config ) ;
	std::cout << "heatup time: " << omp_get_wtime() - time_start << std::endl ;
	energy_compare = ( energy_compare > energy_compare_min ) ? energy_compare : energy_compare_min ;
	T temperature = std::abs( energy_compare - energy_best ) ;

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
	time_start = omp_get_wtime() ;
	time_prev  = time_start ;

	// begin annealing
	#pragma omp parallel
	for( std::size_t iter = 0 ; iter < local_num_iter_max ; ++iter ) {
		
		// test a new configuration, force-driven or random
		move_points_force( step_size ) ;
		//move_points_rand( step_size * 1E-5 ) ;
		update_springs() ;
		update_forces() ; // also computes sum_net_force_magnitude

		#pragma omp single
		{
			// evaluate objective function at new configuration
			energy = get_objective() ;
			// if configuration is too extreme, reset network to best
			reboot = test_reboot() ;
			if( reboot ) {
				std::cout << "REBOOT" << std::endl ;
				points = points_best ;
				points_prev = points_best ;
				energy_prev = energy_best ;
				reboot = false ;
			} else {

				//
				relative_change_energy = std::abs( energy - energy_prev ) / energy_compare ;
				if( relative_change_energy < local_tolerance_change_objective ) {
					++num_small_change ;
				}

				// accept or reject new state based on change in energy
				// accepting decreased energy is guaranteed
				// accepting increased energy is more likely at high temperatures
				if( energy < energy_best ) {
					mean_energy_change += (energy-energy_prev) ;
					step_size /= step_size_reduction ;
					points_best = points ;
					points_prev = points ;
					energy_best = energy ;
					energy_prev = energy ;
					num_consecutive_reject = 0 ;
					num_since_last_best = 0 ;
				} else if( (energy==energy_best) && (sum_net_force_magnitude<sum_net_force_magnitude_prev) ) {
					// why is the final sum_F_net higher than the one just before it????
					points_best = points ;
					points_prev = points ;
					energy_prev = energy ;
					num_consecutive_reject = 0 ;
					num_since_last_best = 0 ;
				} else if( accept_new_points( energy-energy_prev , temperature ) ) {
					mean_energy_change += (energy-energy_prev) ;
					points_prev = points ;
					energy_prev = energy ;
					num_consecutive_reject = 0 ;
					++num_since_last_best ;
				} else {
					step_size *= step_size_reduction ;
					points = points_prev ;
					++num_consecutive_reject ;
					++num_since_last_best ;
				}

				if( (local_num_iter_print>0) && (iter%local_num_iter_print==0) ) {
					mean_energy_change /= static_cast<T>(local_num_iter_print) ;
					time_curr = omp_get_wtime() ;
					std::cout
						<< std::setprecision(3)
						<< std::scientific
						<< "  " << std::setw(10) << energy_best
						<< "  " << std::setw(10) << energy
						<< "  " << std::setw(10) << energy_prev
						<< "  " << std::setw(10) << mean_energy_change
						<< "  " << std::setw(10) << step_size
						<< "  " << std::setw(10) << sum_net_force_magnitude
						<< "  " << std::setw(10) << time_curr - time_prev
						<< std::endl ;
					time_prev = time_curr ;
					mean_energy_change = static_cast<T>(0.0) ;
				}

				if( (local_num_iter_save>0) && (iter%local_num_iter_save==0) ) {
					std::string dir_output_iter_curr = dir_output_iter + "iter_" + std::to_string(iter) + FILESEP ;
					make_dir( dir_output_iter_curr ) ;
					save_network_binary(
						(dir_output_iter_curr+"network_nodes.dat").c_str()   ,
						(dir_output_iter_curr+"network_springs.dat").c_str() ) ;
				}

				// if energy changes are small for multiple iterations, then
				// assume the current state is near local minimum, and
				// reduce the temperature
				if(    ( num_small_change >= num_small_change_max )
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
				break_sum_net_force_magnitude = sum_net_force_magnitude < local_tolerance_sum_net_force ;
				break_consecutive_reject = num_consecutive_reject >= num_consecutive_reject_max ;
				break_isnan_energy = std::isnan( energy ) ;
			}
		}
		if( reboot ) {
			continue ;
		} else if( break_temperature_reductions ) {
			break ;
		} else if( break_sum_net_force_magnitude ) {
			break ;
		} else if( break_consecutive_reject ) {
			break ;
		} else if( break_isnan_energy ) {
			break ;
		}
		/*
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
			points_prev = points ;
			energy_prev = energy ;
		}
		//*/
	}
	std::cout
		<< std::setprecision(3)
		<< std::scientific
		<< "  " << std::setw(10) << energy_best
		<< "  " << std::setw(10) << energy
		<< "  " << std::setw(10) << energy_prev
		<< "  " << std::setw(10) << energy-energy_prev
		<< "  " << std::setw(10) << step_size
		<< "  " << std::setw(10) << sum_net_force_magnitude
		<< std::endl ;

	std::cout << energy - energy_best << std::endl ;

	points = points_best ;
	energy = energy_best ;
	update_springs() ;
	update_forces() ; // also computes sum_net_force_magnitude
	energy = get_objective() ;
	std::cout
		<< std::setprecision(3)
		<< std::scientific
		<< "  " << std::setw(10) << energy_best
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
// ███████  █████  ██    ██ ███████          ██████  ██    ██ ████████ ██████  ██    ██ ████████ 
// ██      ██   ██ ██    ██ ██              ██    ██ ██    ██    ██    ██   ██ ██    ██    ██    
// ███████ ███████ ██    ██ █████           ██    ██ ██    ██    ██    ██████  ██    ██    ██    
//      ██ ██   ██  ██  ██  ██              ██    ██ ██    ██    ██    ██      ██    ██    ██    
// ███████ ██   ██   ████   ███████ ███████  ██████   ██████     ██    ██       ██████     ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::save_output( void )
{
	save_network_binary( file_output_nodes.c_str() , file_output_springs.c_str() ) ;
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

template class spmat<float > ;
template class spmat<double> ;

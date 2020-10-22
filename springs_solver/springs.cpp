//
//  springs.cpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#include "springs.hpp"
#include "vectors_nd.hpp"
#include "klein_summer.hpp"
#include "spmat.hpp"

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

#ifdef _WIN32
	#define FILESEP '\\'
#else
	#define FILESEP '/'
#endif

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


////////////////////////////////////////////////////////////
//////////////////////////  POINT //////////////////////////
////////////////////////////////////////////////////////////

template< class T , std::size_t N >
void Point<T,N>::move( const Vector<T,N> & displacement )
{
	if( !this->fixed_all_dim ) {
		this->position += displacement.zero_where( this->fixed_dim ) ;
	}
}

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
	T force_sign = static_cast<T>(0) ;
	T energy = static_cast<T>(0) ;
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
													  force_magnitude ,
													  energy ) ;
				break ;
			case Spring<T,N>::ForceLengthRelationship::exponential :
				Spring<T,N>::spring_force_exponential( delta_length ,
													   * num_force_length_parameters ,
													   * force_length_parameters ,
													   force_magnitude ,
													   energy ) ;
				break ;
			case Spring<T,N>::ForceLengthRelationship::powerlaw :
				Spring<T,N>::spring_force_powerlaw( delta_length ,
												    * num_force_length_parameters ,
												    * force_length_parameters ,
												    force_magnitude ,
												    energy ) ;
				break ;
			default:
				force_magnitude = static_cast<T>(0) ;
				energy = static_cast<T>(0) ;
		}
		effective_spring_constant = force_magnitude / delta_length ;
		force_magnitude *= force_sign ;
	}
	force = delta_position ;
	force *= (force_magnitude/length) ;
	return energy ;
}

///______________  spring_force_polynomial ______________///
template< class T , std::size_t N >
void Spring<T,N>::spring_force_polynomial( const T & delta_length ,
										   const std::size_t & num_force_length_parameters ,
										   const std::vector<T> & force_length_parameters ,
										   T & force_magnitude ,
										   T & energy )
{
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T force_component ;
	T delta_length_power_i = delta_length ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; ++i ) {
		force_component = force_length_parameters[i] * delta_length_power_i ;
		force_magnitude += force_component ;
		energy += force_component * delta_length / static_cast<T>(i+1) ;
		delta_length_power_i *= delta_length ;
	}
	return ;
}

///______________ spring_force_exponential ______________///
template< class T , std::size_t N >
void Spring<T,N>::spring_force_exponential( const T & delta_length ,
										    const std::size_t & num_force_length_parameters ,
										    const std::vector<T> & force_length_parameters ,
										    T & force_magnitude ,
										    T & energy )
{
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T force_component ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; i+=2 ) {
		force_component = force_length_parameters[i] * std::expm1( delta_length * force_length_parameters[i+1] ) ;
		force_magnitude += force_component ;
		energy += (force_component/force_length_parameters[i+1]) + (delta_length*force_length_parameters[i]) ;
	}
	return ;
}

///_______________  spring_force_powerlaw _______________///
template< class T , std::size_t N >
void Spring<T,N>::spring_force_powerlaw( const T & delta_length ,
										 const std::size_t & num_force_length_parameters ,
										 const std::vector<T> & force_length_parameters ,
										 T & force_magnitude ,
										 T & energy )
{
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T force_component ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; i+=2 ) {
		force_component = force_length_parameters[i] * std::pow( delta_length , force_length_parameters[i+1] ) ;
		force_magnitude += force_component ;
		energy += force_component * delta_length / (force_length_parameters[i+1]+static_cast<T>(1)) ;
	}
	return ;
}

///_______________  spring_stiffness_rest _______________///
template< class T , std::size_t N >
T Spring<T,N>::spring_stiffness_rest( void )
{
	std::size_t * num_force_length_parameters ;
	std::vector<T> * force_length_parameters ;
	T effective_spring_constant_rest = static_cast<T>(0) ;
	if( num_force_length_parameters_tension > 0 ) {
		num_force_length_parameters = & num_force_length_parameters_tension ;
		force_length_parameters     = & force_length_parameters_tension ;
	} else if( num_force_length_parameters_compression > 0 ) {
		num_force_length_parameters = & num_force_length_parameters_compression ;
		force_length_parameters     = & force_length_parameters_compression ;
	} else {
		return static_cast<T>(0) ;
	}
	switch( force_length_type_tension ) {
		case Spring<T,N>::ForceLengthRelationship::polynomial :
			effective_spring_constant_rest = Spring<T,N>::spring_stiffness_rest_polynomial( num_force_length_parameters_tension ,
																							force_length_parameters_tension ) ;
			break ;
		case Spring<T,N>::ForceLengthRelationship::exponential :
			effective_spring_constant_rest = Spring<T,N>::spring_stiffness_rest_exponential( num_force_length_parameters_tension ,
																							 force_length_parameters_tension ) ;
			break ;
		case Spring<T,N>::ForceLengthRelationship::powerlaw :
			effective_spring_constant_rest = Spring<T,N>::spring_stiffness_rest_powerlaw( num_force_length_parameters_tension ,
																						  force_length_parameters_tension ) ;
			break ;
		default:
			effective_spring_constant_rest = static_cast<T>(0) ;
	}
	return effective_spring_constant_rest ;
}

///__________ spring_stiffness_rest_polynomial __________///
template< class T , std::size_t N >
T Spring<T,N>::spring_stiffness_rest_polynomial( const std::size_t & num_force_length_parameters ,
												 const std::vector<T> & force_length_parameters )
{
	return force_length_parameters[0] ;
}

///_________  spring_stiffness_rest_exponential _________///
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

///___________ spring_stiffness_rest_powerlaw ___________///
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

///_____________ spring_get_scale_stiffness _____________///
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

///______________ spring_rescale_stiffness ______________///
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

/// DECLARE STATIC MEMBER VARIABLES ///
// none to declare

////////////////////////////////////////////////////////////
////////////////////// SPRING NETWORK //////////////////////
////////////////////////////////////////////////////////////

///_____________  create_spring_network_obj _____________///
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

///_______________________  setup _______________________///
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
		read_data_uint32  = new std::uint32_t [2] ;
		read_data_uint8   = new std::uint8_t  [4] ;
		for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
			std::fread( read_data_uint32 , sizeof(std::uint32_t) , 2 , file_ptr ) ;
			std::fread( read_data_uint8  , sizeof(std::uint8_t)  , 4 , file_ptr ) ;
			NFLPT = read_data_uint8[2] ;
			NFLPC = read_data_uint8[3] ;
			read_data_T = new T [1+NFLPT+NFLPC] ;
			std::fread( read_data_T , sizeof(T) , 1+NFLPT+NFLPC , file_ptr ) ;
			s->start = &points[ read_data_uint32[0] ] ;
			s->end   = &points[ read_data_uint32[1] ] ;
			s->force_length_type_tension     = static_cast<typename Spring<T,N>::ForceLengthRelationship>(read_data_uint8[0]) ;
			s->force_length_type_compression = static_cast<typename Spring<T,N>::ForceLengthRelationship>(read_data_uint8[1]) ;
			s->num_force_length_parameters_tension     = NFLPT ;
			s->num_force_length_parameters_compression = NFLPC ;
			s->rest_length = read_data_T[0] ;
			s->force_length_parameters_tension     = std::vector<T>( &read_data_T[1      ] , &read_data_T[1+NFLPT      ] ) ;
			s->force_length_parameters_compression = std::vector<T>( &read_data_T[1+NFLPT] , &read_data_T[1+NFLPT+NFLPC] ) ;
			delete [] read_data_T ;
		}
		delete [] read_data_uint32 ;
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
		write_data_uint32 = new std::uint32_t [2] ;
		write_data_uint8  = new std::uint8_t  [4] ;
		for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
			NFLPT = s->num_force_length_parameters_tension ;
			NFLPC = s->num_force_length_parameters_compression ;
			write_data_T = new T [1+NFLPT+NFLPC] ;
			write_data_uint32[0] = static_cast<std::uint32_t>( s->start - &points[0] ) ;
			write_data_uint32[1] = static_cast<std::uint32_t>( s->end   - &points[0] ) ;
			write_data_uint8[0] = static_cast<std::uint8_t>( s->force_length_type_tension     ) ;
			write_data_uint8[1] = static_cast<std::uint8_t>( s->force_length_type_compression ) ;
			write_data_uint8[2] = static_cast<std::uint8_t>( NFLPT ) ;
			write_data_uint8[3] = static_cast<std::uint8_t>( NFLPC ) ;
			write_data_T[0] = s->rest_length ;
			std::copy( s->force_length_parameters_tension.begin()     , s->force_length_parameters_tension.end()     , &write_data_T[1      ] ) ;
			std::copy( s->force_length_parameters_compression.begin() , s->force_length_parameters_compression.end() , &write_data_T[1+NFLPT] ) ;
			std::fwrite( write_data_uint32 , sizeof(std::uint32_t) , 2             , file_ptr ) ;
			std::fwrite( write_data_uint8  , sizeof(std::uint8_t)  , 4             , file_ptr ) ;
			std::fwrite( write_data_T      , sizeof(T)             , 1+NFLPT+NFLPC , file_ptr ) ;
			delete [] write_data_T ;
		}
		delete [] write_data_uint32 ;
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
		points[p].not_referenced = true ;
	}
	std::size_t p_start ;
	std::size_t p_end ;
	for( std::size_t s = 0 ; s < num_springs ; ++s ) {
		p_start = springs[s].start - &points[0] ;
		p_end   = springs[s].end   - &points[0] ;
		nodes[p_start].links.push_back( (Link){ static_cast<std::size_t>(0) , &points[p_end  ] , &springs[s] , static_cast<T>(+1) } ) ;
		nodes[p_end  ].links.push_back( (Link){ static_cast<std::size_t>(0) , &points[p_start] , &springs[s] , static_cast<T>(-1) } ) ;
		points[p_start].not_referenced = false ;
		points[p_end  ].not_referenced = false ;
	}
	nodes.erase( std::remove_if( nodes.begin() , nodes.end() , []( Node n ){ return n.point->fixed_all_dim  ; } ) , nodes.end() ) ;
	nodes.erase( std::remove_if( nodes.begin() , nodes.end() , []( Node n ){ return n.point->not_referenced ; } ) , nodes.end() ) ;
	for( std::size_t i = 0 ; i < nodes.size() ; ++i ) {
		nodes[i].node_index = i ;
	}
	iterNode n ;
	iterNode n_end ;
	iterLink l ;
	iterLink l_end ;
	std::size_t i ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			for( i = 0 ; i < nodes.size() ; ++i ) {
				if( nodes[i].point == l->point ) {
					break ;
				}
			}
			l->node_index = i ;
		}
	}
	return ;
}

///_______________________  solve _______________________///
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

///_________________  get_scale_network _________________///
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

///__________________  rescale_network __________________///
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

///_______________ find_max_spring_length _______________///
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

///____________________ total_energy ____________________///
template< class T , std::size_t N >
T SpringNetwork<T,N>::total_energy( void )
{
	// using klein summation, modification of kahan summation,
	// to account for floating point precision issues with sums
	// of many small numbers, or small + large numbers
	KleinSummer<T> ksum ;

	// internal energy due to spring tension/compression
	ksum.reset() ;
	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
		ksum.add( s->spring_energy() ) ;
	}
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
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		p->net_force = static_cast<T>(0.0) ;
	}
	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
		s->start->net_force += s->get_force() ;
		s->end->net_force   -= s->get_force() ;
	}
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		p->net_force += p->force_applied ;
		p->net_force_magnitude = p->net_force.norm() ;
	}
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

	//energy += ( sum_net_force_magnitude * scale_length ) ;

	if( objective=="energy"   ) { return energy                  ; } else
	if( objective=="sumforce" ) { return sum_net_force_magnitude ; } else
	if( objective=="maxforce" ) { return max_net_force_magnitude ; }
	else { return static_cast<T>(0.0) ; }
}

///__________________ compute_gradient __________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_gradient( void )
{
	// evalute net force at each node (gradient of total spring network energy)
	std::size_t num_node = nodes.size() ;
	neg_gradient.resize( N*num_node ) ;
	total_energy() ;
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

///_____________  compute_hessian_numerical _____________///
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
			if( l->node_index < nodes.size() ) {
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

///_____________ compute_hessian_analytical _____________///
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

///___________  compute_newton_step_direction ___________///
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
	if( dot_product <= 0.0 ) {
		step_direction = neg_gradient ;
	}
	//*/

	return ;
}

///_________________ move_points_newton _________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_newton( const T & step_size )
{
	// move small displacement towards equilibrating position by 2nd order newton method
	iterNode n ;
	iterNode n_end ;
	std::size_t d ;
	Vector<T,N> small_displacement ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( d = 0 ; d < N ; ++d ) {
			small_displacement[d] = step_direction[ (n->node_index)*N + d ] ;
		}
		small_displacement *= step_size ;
		n->point->move( small_displacement ) ;
	}
	return ;
}

///_________________  move_points_force _________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_force( const T & step_size )
{
	// apply net force & move small displacement towards equilibrating position
	iterNode n ;
	iterNode n_end ;
	Vector<T,N> small_displacement ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		small_displacement = n->point->net_force ;
		small_displacement *= step_size ;
		n->point->move( small_displacement ) ;
	}
	return ;
}

///__________________ move_points_rand __________________///
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

///__________________ minimize_energy __________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::minimize_energy( void )
{
	// copy current state to initial, previous, and best states
	points_init = points ;
	points_prev = points ;
	T energy = total_energy() ;
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

	//
	std::cout
		<< "  " << std::setw(10) << "E_curr"
		<< "  " << std::setw(10) << "E_prev"
		<< "  " << std::setw(10) << "E_change"
		<< "  " << std::setw(10) << "stepSize"
		<< "  " << std::setw(10) << "sum_F_net"
		<< std::endl ;

	//
	for( std::size_t iter = 0 ; iter < local_num_iter_max ; ++iter ) {

		// perturb the system with some randomness
		if( is_shaking ) {
			if( iter%num_iter_shake == 0 ) {
				move_points_rand( shake_step_size ) ;
				shake_step_size *= shake_step_size_reduction ;
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
			energy = total_energy() ; // also computes sum_net_force_magnitude
		} while( (energy>energy_prev) && (step_size>step_size_min) ) ;
		/*
		while( energy == energy_prev ) {
			step_size /= step_size_reduction ;
			points = points_prev ;
			move_points_force( step_size ) ;
			energy = total_energy() ; // also computes sum_net_force_magnitude
		}
		//*/

		//
		change_energy = (energy_prev-energy) / (energy_init-energy_prev) ;
		if( energy == energy_prev ) {
			++num_iter_zero_change ;
		}

		if( (local_num_iter_print>0) && (iter%local_num_iter_print==0) ) {
			std::cout
				<< std::setprecision(3)
				<< std::scientific
				<< "  " << std::setw(10) << energy
				<< "  " << std::setw(10) << energy_prev
				<< "  " << std::setw(10) << energy-energy_prev
				<< "  " << std::setw(10) << step_size
				<< "  " << std::setw(10) << sum_net_force_magnitude
				<< std::endl ;
		}

		if( (local_num_iter_save>0) && (iter%local_num_iter_save==0) ) {
			std::string dir_output_iter_curr = dir_output_iter + "iter_" + std::to_string(iter) + FILESEP ;
			make_dir( dir_output_iter_curr ) ;
			save_network_binary(
				(dir_output_iter_curr+"network_nodes.dat").c_str()   ,
				(dir_output_iter_curr+"network_springs.dat").c_str() ) ;
		}

		// stopping condition
		if(    (sum_net_force_magnitude<local_tolerance_sum_net_force)
			|| (change_energy<local_tolerance_change_objective)
			|| (step_size<=step_size_min)
			|| (num_iter_zero_change>num_iter_zero_change_max) )
		{
			break ;
		}
		if( std::isnan( energy ) ) {
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
	return ;
}

///______________  minimize_energy_newton ______________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::minimize_energy_newton( void )
{
	// copy current state to initial, previous, and best states
	points_init = points ;
	points_prev = points ;
	T energy = total_energy() ;
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

	//
	std::cout
		<< "  " << std::setw(10) << "E_curr"
		<< "  " << std::setw(10) << "E_prev"
		<< "  " << std::setw(10) << "E_change"
		<< "  " << std::setw(10) << "stepSize"
		<< "  " << std::setw(10) << "sum_F_net"
		<< std::endl ;

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
			energy = total_energy() ; // also computes sum_net_force_magnitude
		} while( (energy>=energy_prev) && (step_size>step_size_min) ) ;

		//
		change_energy = (energy_prev-energy) / (energy_init-energy_prev) ;
		if( energy == energy_prev ) {
			++num_iter_zero_change ;
		}

		if( (local_num_iter_print>0) && (iter%local_num_iter_print==0) ) {
			std::cout
				<< std::setprecision(3)
				<< std::scientific
				<< "  " << std::setw(10) << energy
				<< "  " << std::setw(10) << energy_prev
				<< "  " << std::setw(10) << energy-energy_prev
				<< "  " << std::setw(10) << step_size
				<< "  " << std::setw(10) << sum_net_force_magnitude
				<< std::endl ;
		}

		if( (local_num_iter_save>0) && (iter%local_num_iter_save==0) ) {
			std::string dir_output_iter_curr = dir_output_iter + "iter_" + std::to_string(iter) + FILESEP ;
			make_dir( dir_output_iter_curr ) ;
			save_network_binary(
				(dir_output_iter_curr+"network_nodes.dat").c_str()   ,
				(dir_output_iter_curr+"network_springs.dat").c_str() ) ;
		}

		// stopping condition
		if(    (sum_net_force_magnitude<local_tolerance_sum_net_force)
			|| (change_energy<local_tolerance_change_objective)
			|| (step_size<=step_size_min)
			|| (num_iter_zero_change>num_iter_zero_change_max) )
		{
			break ;
		}
		if( std::isnan( energy ) ) {
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
	//T energy_init = energy ;
	T energy_prev = energy ;
	T energy_best = energy ;
	
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
	find_max_spring_length() ;
	
	// choose starting temperature
	// equivalent to 1/e probability to accept energy difference
	T heatup_amplitude = 0.01 ;
	std::size_t heatup_num_config = 1000 ;
	energy_compare = heat_up( heatup_amplitude , heatup_num_config ) ;
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
		<< std::endl ;

	// begin annealing
	for( std::size_t iter = 0 ; iter < local_num_iter_max ; ++iter ) {
		
		// test a new configuration, force-driven or random
		move_points_force( step_size ) ;
		//move_points_rand( step_size * 1E-5 ) ;
		energy = total_energy() ; // also computes sum_net_force_magnitude

		// if configuration is too extreme, reset network to best
		reboot = test_reboot() ;
		if( reboot ) {
			std::cout << "REBOOT" << std::endl ;
			points = points_best ;
			points_prev = points_best ;
			energy_prev = energy_best ;
			reboot = false ;
			continue ;
		}

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
		}
		else if( accept_new_points( energy-energy_prev , temperature ) ) {
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
			std::cout
				<< std::setprecision(3)
				<< std::scientific
				<< "  " << std::setw(10) << energy_best
				<< "  " << std::setw(10) << energy
				<< "  " << std::setw(10) << energy_prev
				<< "  " << std::setw(10) << mean_energy_change
				<< "  " << std::setw(10) << step_size
				<< "  " << std::setw(10) << sum_net_force_magnitude
				<< std::endl ;
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
				break ;
			}
			num_small_change = 0 ;
			num_since_last_best = 0 ;
		}
		//
		if( sum_net_force_magnitude < local_tolerance_sum_net_force ) {
			break ;
		}
		if( num_consecutive_reject >= num_consecutive_reject_max ) {
			break ;
		}
		if( std::isnan( energy ) ) {
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
	points = points_best ;
	energy = energy_best ;
	return ;
}

///____________________  save_output ____________________///
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

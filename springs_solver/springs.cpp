//
//  springs.cpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#include "springs.hpp"
#include "vectors_nd.hpp"
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
				if( arg=="num_points"                ) { num_points                = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_springs"               ) { num_springs               = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="precision"                 ) { precision                 =                                    val   ; } else
				if( arg=="num_dimensions"            ) { num_dimensions            = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_stiffness_tension"     ) { num_stiffness_tension     = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_stiffness_compression" ) { num_stiffness_compression = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="algorithm"                 ) { algorithm                 =                                    val   ; } else
				if( arg=="num_iter_save"             ) { num_iter_save             = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_print"            ) { num_iter_print            = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_max"              ) { num_iter_max              = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="tolerance_change_energy"   ) { tolerance_change_energy   =                          std::stod(val)  ; } else
				if( arg=="tolerance_sum_net_force"   ) { tolerance_sum_net_force   =                          std::stod(val)  ; }
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
			file << "num_points"                << ' ' << num_points                << '\n'
				 << "num_springs"               << ' ' << num_springs               << '\n'
				 << "precision"                 << ' ' << precision                 << '\n'
				 << "num_dimensions"            << ' ' << num_dimensions            << '\n'
				 << "num_stiffness_tension"     << ' ' << num_stiffness_tension     << '\n'
				 << "num_stiffness_compression" << ' ' << num_stiffness_compression << '\n'
				 << "algorithm"                 << ' ' << algorithm                 << '\n'
				 << "num_iter_save"             << ' ' << num_iter_save             << '\n'
				 << "num_iter_print"            << ' ' << num_iter_print            << '\n'
				 << "num_iter_max"              << ' ' << num_iter_max              << '\n'
				 << "tolerance_change_energy"   << ' ' << tolerance_change_energy   << '\n'
				 << "tolerance_sum_net_force"   << ' ' << tolerance_sum_net_force ;
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
	//
	algorithm               = network_parameters.algorithm ;
	num_iter_save           = network_parameters.num_iter_save ;
	num_iter_print          = network_parameters.num_iter_print ;
	num_iter_max            = network_parameters.num_iter_max ;
	tolerance_change_energy = network_parameters.tolerance_change_energy ;
	tolerance_sum_net_force = network_parameters.tolerance_sum_net_force ;
	//
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
	nodes.erase( std::remove_if( nodes.begin() , nodes.end() , []( Node n ){ return n.point->fixed          ; } ) , nodes.end() ) ;
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
		scale_stiffness += s->stiffness_tension[0] ;
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
	switch( direction ) {
		case +1:
			// set to normalize scale from original scale
			for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
				p->position -= scale_position_min ;
				p->position /= scale_length ;
				p->force_applied /= scale_force ;
			}
			for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
				s->rest_length /= scale_length ;
				for( std::size_t i = 0 ; i < Spring<T,N>::num_stiffness_tension ; ++i ) {
					s->stiffness_tension[i] /= scale_stiffness ;
				}
				for( std::size_t i = 0 ; i < Spring<T,N>::num_stiffness_compression ; ++i ) {
					s->stiffness_compression[i] /= scale_stiffness ;
				}
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
				s->rest_length *= scale_length ;
				for( std::size_t i = 0 ; i < Spring<T,N>::num_stiffness_tension ; ++i ) {
					s->stiffness_tension[i] *= scale_stiffness ;
				}
				for( std::size_t i = 0 ; i < Spring<T,N>::num_stiffness_compression ; ++i ) {
					s->stiffness_compression[i] *= scale_stiffness ;
				}
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
	energy += ( sum_net_force_magnitude * scale_length ) ;

	return energy ;
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

///__________________  compute_hessian __________________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_hessian( void )
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
		// change in net force at n'th node w.r.t. small changes in linked node positions
		// only consider moveable nodes (i.e., not fixed points)
		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
			if( l->node_index < nodes.size() ) {
				for( d_p = 0 ; d_p < N ; ++d_p ) {
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
	hessian = typename spmat<T>::spmat( N*num_node , N*num_node ) ;
	hessian.set_val( sp_row , sp_col , sp_val ) ;
	return ;
}

///___________  compute_newton_step_direction ___________///
template< class T , std::size_t N >
void SpringNetwork<T,N>::compute_newton_step_direction( void )
{
	// prepare gradient and hessian
	compute_gradient() ;
	compute_hessian() ;
	step_direction = neg_gradient ;

	// if any zero on the diagonal, remove row and column
	// set all row and col values 0, set diagonal value 1, set source vector 0
	T small_number = static_cast<T>(1000) * std::numeric_limits<T>::epsilon() ;
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
		n->point->position += small_displacement ;
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
		n->point->position += small_displacement ;
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
		n->point->position += small_displacement ;
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
	T step_size_min = 1E-16 ;
	T change_energy ;

	//
	bool is_shaking = false ;
	std::size_t num_iter_shake = 100 ;
	T shake_step_size = 1E-2 ;
	T shake_step_size_reduction = 0.999 ;

	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 500000 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :   5000 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :    200 ;
	T local_tolerance_sum_net_force = ( tolerance_sum_net_force > 0.0 ) ? tolerance_sum_net_force : 1E-12 ;
	T local_tolerance_change_energy = ( tolerance_change_energy > 0.0 ) ? tolerance_change_energy : 1E-12 ;

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
		} while( (energy>=energy_prev) && (step_size>step_size_min) ) ;

		//
		change_energy = (energy_prev-energy) / (energy_init-energy_prev) ;
		if( energy == energy_prev ) {
			++num_iter_zero_change ;
		}

		if( iter%local_num_iter_print == 0 ) {
			std::cout
				<< std::setprecision(3)
				<< std::scientific
				<< "  " << std::setw(10) << energy
				<< "  " << std::setw(10) << energy_prev
				<< "  " << std::setw(10) << energy-energy_prev
				<< "  " << std::setw(10) << step_size
				<< std::endl ;
		}

		// stopping condition
		if(    (sum_net_force_magnitude<local_tolerance_sum_net_force)
			|| (change_energy<local_tolerance_change_energy)
			|| (step_size<=step_size_min)
			|| (num_iter_zero_change>num_iter_zero_change_max) )
		{
			break ;
		}
	}
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
	T step_size = 1E-6 ;
	T step_size_reduction = 0.9 ; // in range (0,1)
	T step_size_min = 1E-16 ;
	T change_energy ;

	// user-definable parameters with default values specific to this algorithm
	std::size_t local_num_iter_max   = ( num_iter_max   > 0 ) ? num_iter_max   : 500 ;
	std::size_t local_num_iter_save  = ( num_iter_save  > 0 ) ? num_iter_save  :  10 ;
	std::size_t local_num_iter_print = ( num_iter_print > 0 ) ? num_iter_print :   1 ;
	T local_tolerance_sum_net_force = ( tolerance_sum_net_force > 0.0 ) ? tolerance_sum_net_force : 1E-12 ;
	T local_tolerance_change_energy = ( tolerance_change_energy > 0.0 ) ? tolerance_change_energy : 1E-12 ;

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
			points = points_prev ;
			move_points_newton( step_size ) ;
			energy = total_energy() ; // also computes sum_net_force_magnitude
		} while( (energy>=energy_prev) && (step_size>step_size_min) ) ;

		//
		change_energy = (energy_prev-energy) / (energy_init-energy_prev) ;
		if( energy == energy_prev ) {
			++num_iter_zero_change ;
		}

		if( iter%local_num_iter_print == 0 ) {
			std::cout
				<< std::setprecision(3)
				<< std::scientific
				<< "  " << std::setw(10) << energy
				<< "  " << std::setw(10) << energy_prev
				<< "  " << std::setw(10) << energy-energy_prev
				<< "  " << std::setw(10) << step_size
				<< std::endl ;
		}

		// stopping condition
		if(    (sum_net_force_magnitude<local_tolerance_sum_net_force)
			|| (change_energy<local_tolerance_change_energy)
			|| (step_size<=step_size_min)
			|| (num_iter_zero_change>num_iter_zero_change_max) )
		{
			break ;
		}
	}
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
	// TODO // get these parameters from NetworkParameters? (default/file)
	T step_size = static_cast<T>(0.02) ;
	T temperature_reduction = static_cast<T>(0.99) ;
	T energy_compare ;
	T energy_compare_min = static_cast<T>(1E-12) ;
	T relative_change_energy ;
	T relative_change_energy_tol = static_cast<T>(1E-6) ;
	std::size_t num_iter_max       = 500000 ;
	std::size_t num_iter_heatup    = num_iter_max /   5 ;
	std::size_t num_iter_10percent = num_iter_max /  10 ;
	std::size_t num_iter_01percent = num_iter_max / 100 ;
	std::size_t num_consecutive_reject     =  0 ;
	std::size_t num_consecutive_reject_max = 20 ;
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

		//	print status
		if( iter%num_iter_01percent == 0 ) {
			if( iter > 0 ) {
				std::cout << "." << std::flush ;
				if( iter%num_iter_10percent == 0 ) {
					std::cout << 100*iter/num_iter_max << "%\n" << std::flush ;
				}
			} else {
				std::cout << '\n' << std::flush ;
			}
		}
		
		// randomize the order of nodes?
		// no need for random ordering if force/energy does not depend on order
		///std::random_shuffle( nodes.begin() , nodes.end() ) ;
		
		// test a new configuration, force-driven or random
		move_points_force( step_size ) ;
		//move_points_rand( step_size * 1E-5 ) ;
		energy = total_energy() ;
		
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
		
		// if energy changes are small for multiple iterations, then
		// assume the current state is near local minimum, and
		// reduce the temperature
		relative_change_energy = std::fabs( energy - energy_prev ) / energy_compare ;
		//**
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

		//*/
		if( num_consecutive_reject >= num_consecutive_reject_max ) {
			num_consecutive_reject = 0 ;
			step_size *= temperature_reduction ; // TEST REMOVE
		}
		//**/
		
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
			num_consecutive_reject = 0 ;
		}
		else if( accept_new_points( energy-energy_prev , temperature ) ) {
			points_prev = points ;
			energy_prev = energy ;
			num_consecutive_reject = 0 ;
		} else {
			points = points_prev ;
			++num_consecutive_reject ;
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

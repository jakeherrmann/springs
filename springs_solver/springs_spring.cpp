#include "springs.hpp"
#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ██████  ███████  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████         ██████   █████  ██████   █████  ███    ███ ███████ ████████ ███████ ██████  ███████ 
//  ██   ██ ██   ██ ██      ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██              ██   ██ ██   ██ ██   ██ ██   ██ ████  ████ ██         ██    ██      ██   ██ ██      
//  ██████  ██████  █████   ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████           ██████  ███████ ██████  ███████ ██ ████ ██ █████      ██    █████   ██████  ███████ 
//  ██      ██   ██ ██      ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██              ██      ██   ██ ██   ██ ██   ██ ██  ██  ██ ██         ██    ██      ██   ██      ██ 
//  ██      ██   ██ ███████  ██████  ██████  ██      ██ ██       ██████     ██    ███████ ███████ ██      ██   ██ ██   ██ ██   ██ ██      ██ ███████    ██    ███████ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::precompute_parameters( void )
{
	// tension
	switch( force_length_type_tension ) {
		//
		case Spring<T,N>::ForceLengthRelationship::polynomial :
			energy_length_parameters_tension.resize( num_force_length_parameters_tension ) ;
			stiffness_length_parameters_tension.resize( num_force_length_parameters_tension ) ;
			for( std::size_t i = 0 ; i < num_force_length_parameters_tension ; ++i ) {
				energy_length_parameters_tension   [i] = force_length_parameters_tension[i] / static_cast<T>(i+2) ;
				stiffness_length_parameters_tension[i] = force_length_parameters_tension[i] / static_cast<T>(i+1) ;
			}
			break ;
		//
		case Spring<T,N>::ForceLengthRelationship::exponential :
			energy_length_parameters_tension.resize( num_force_length_parameters_tension ) ;
			for( std::size_t i = 0 ; i < num_force_length_parameters_tension ; i+=2 ) {
				energy_length_parameters_tension[i] = static_cast<T>(1.0) / force_length_parameters_tension[i+1] ;
			}
			break ;
		//	
		case Spring<T,N>::ForceLengthRelationship::powerlaw :
			energy_length_parameters_tension.resize( num_force_length_parameters_tension ) ;
			stiffness_length_parameters_tension.resize( num_force_length_parameters_tension ) ;
			for( std::size_t i = 0 ; i < num_force_length_parameters_tension ; i+=2 ) {
				energy_length_parameters_tension   [i] = force_length_parameters_tension[i] / (force_length_parameters_tension[i+1]+static_cast<T>(1.0)) ;
				stiffness_length_parameters_tension[i] = force_length_parameters_tension[i] /  force_length_parameters_tension[i+1] ;
			}
			break ;
		default: ;
	}
	// compression
	switch( force_length_type_compression ) {
		//
		case Spring<T,N>::ForceLengthRelationship::polynomial :
			energy_length_parameters_compression.resize( num_force_length_parameters_compression ) ;
			stiffness_length_parameters_compression.resize( num_force_length_parameters_compression ) ;
			for( std::size_t i = 0 ; i < num_force_length_parameters_compression ; ++i ) {
				energy_length_parameters_compression   [i] = force_length_parameters_compression[i] / static_cast<T>(i+2) ;
				stiffness_length_parameters_compression[i] = force_length_parameters_compression[i] / static_cast<T>(i+1) ;
			}
			break ;
		//
		case Spring<T,N>::ForceLengthRelationship::exponential :
			energy_length_parameters_compression.resize( num_force_length_parameters_compression ) ;
			for( std::size_t i = 0 ; i < num_force_length_parameters_compression ; i+=2 ) {
				energy_length_parameters_compression[i] = static_cast<T>(1.0) / force_length_parameters_compression[i+1] ;
			}
			break ;
		//	
		case Spring<T,N>::ForceLengthRelationship::powerlaw :
			energy_length_parameters_compression.resize( num_force_length_parameters_compression ) ;
			stiffness_length_parameters_compression.resize( num_force_length_parameters_compression ) ;
			for( std::size_t i = 0 ; i < num_force_length_parameters_compression ; i+=2 ) {
				energy_length_parameters_compression   [i] = force_length_parameters_compression[i] / (force_length_parameters_compression[i+1]+static_cast<T>(1.0)) ;
				stiffness_length_parameters_compression[i] = force_length_parameters_compression[i] /  force_length_parameters_compression[i+1] ;
			}
			break ;
		default: ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ███    ██ ███████ ██████   ██████  ██    ██ 
//  ██      ████   ██ ██      ██   ██ ██        ██  ██  
//  █████   ██ ██  ██ █████   ██████  ██   ███   ████   
//  ██      ██  ██ ██ ██      ██   ██ ██    ██    ██    
//  ███████ ██   ████ ███████ ██   ██  ██████     ██    
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
	std::vector<T> * energy_length_parameters = NULL ;
	std::vector<T> * stiffness_length_parameters = NULL ;
	if( delta_length == 0.0 ) {
		force_magnitude = static_cast<T>(0) ;
		energy = static_cast<T>(0) ;
		effective_spring_constant = spring_stiffness_rest() ;
	} else {
		if( delta_length > 0 ) {
			force_length_type           = & force_length_type_tension ;
			num_force_length_parameters = & num_force_length_parameters_tension ;
			force_length_parameters     = & force_length_parameters_tension ;
			energy_length_parameters    = & energy_length_parameters_tension ;
			stiffness_length_parameters = & stiffness_length_parameters_tension ;
			force_sign = static_cast<T>(+1) ;
		} else if( delta_length < 0 ) {
			delta_length = std::abs( delta_length ) ;
			force_length_type           = & force_length_type_compression ;
			num_force_length_parameters = & num_force_length_parameters_compression ;
			force_length_parameters     = & force_length_parameters_compression ;
			energy_length_parameters    = & energy_length_parameters_compression ;
			stiffness_length_parameters = & stiffness_length_parameters_compression ;
			force_sign = static_cast<T>(-1) ;
		}
		switch( * force_length_type ) {
			case Spring<T,N>::ForceLengthRelationship::polynomial :
				Spring<T,N>::spring_force_polynomial( delta_length ,
													  * num_force_length_parameters ,
													  * force_length_parameters ,
													  * energy_length_parameters ,
													  * stiffness_length_parameters ,
													  effective_spring_constant ,
													  force_magnitude ,
													  energy ) ;
				break ;
			case Spring<T,N>::ForceLengthRelationship::exponential :
				Spring<T,N>::spring_force_exponential( delta_length ,
													   * num_force_length_parameters ,
													   * force_length_parameters ,
													   * energy_length_parameters ,
													   * stiffness_length_parameters ,
													   effective_spring_constant ,
													   force_magnitude ,
													   energy ) ;
				break ;
			case Spring<T,N>::ForceLengthRelationship::powerlaw :
				Spring<T,N>::spring_force_powerlaw( delta_length ,
												    * num_force_length_parameters ,
												    * force_length_parameters ,
													* energy_length_parameters ,
													* stiffness_length_parameters ,
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
//  ███████  ██████  ██████   ██████ ███████         ██████   ██████  ██      ██    ██ ███    ██  ██████  ███    ███ ██  █████  ██      
//  ██      ██    ██ ██   ██ ██      ██              ██   ██ ██    ██ ██       ██  ██  ████   ██ ██    ██ ████  ████ ██ ██   ██ ██      
//  █████   ██    ██ ██████  ██      █████           ██████  ██    ██ ██        ████   ██ ██  ██ ██    ██ ██ ████ ██ ██ ███████ ██      
//  ██      ██    ██ ██   ██ ██      ██              ██      ██    ██ ██         ██    ██  ██ ██ ██    ██ ██  ██  ██ ██ ██   ██ ██      
//  ██       ██████  ██   ██  ██████ ███████ ███████ ██       ██████  ███████    ██    ██   ████  ██████  ██      ██ ██ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::spring_force_polynomial( const T & delta_length ,
										   const std::size_t & num_force_length_parameters ,
										   const std::vector<T> & force_length_parameters ,
										   const std::vector<T> & energy_length_parameters ,
										   const std::vector<T> & stiffness_length_parameters ,
										   T & effective_spring_constant ,
										   T & force_magnitude ,
										   T & energy )
{
	effective_spring_constant = static_cast<T>(0) ;
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T delta_length_power_im1 = static_cast<T>(1) ;
	T delta_length_power_i = delta_length ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; ++i ) {
		effective_spring_constant += stiffness_length_parameters[i] * delta_length_power_im1 ;
		force_magnitude += force_length_parameters[i] * delta_length_power_i ;
		delta_length_power_im1 = delta_length_power_i ;
		delta_length_power_i *= delta_length ;
		energy += energy_length_parameters[i] * delta_length_power_i ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████  ██████  ██████   ██████ ███████         ███████ ██   ██ ██████   ██████  ███    ██ ███████ ███    ██ ████████ ██  █████  ██      
//  ██      ██    ██ ██   ██ ██      ██              ██       ██ ██  ██   ██ ██    ██ ████   ██ ██      ████   ██    ██    ██ ██   ██ ██      
//  █████   ██    ██ ██████  ██      █████           █████     ███   ██████  ██    ██ ██ ██  ██ █████   ██ ██  ██    ██    ██ ███████ ██      
//  ██      ██    ██ ██   ██ ██      ██              ██       ██ ██  ██      ██    ██ ██  ██ ██ ██      ██  ██ ██    ██    ██ ██   ██ ██      
//  ██       ██████  ██   ██  ██████ ███████ ███████ ███████ ██   ██ ██       ██████  ██   ████ ███████ ██   ████    ██    ██ ██   ██ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::spring_force_exponential( const T & delta_length ,
										    const std::size_t & num_force_length_parameters ,
										    const std::vector<T> & force_length_parameters ,
										    const std::vector<T> & energy_length_parameters ,
										    const std::vector<T> & stiffness_length_parameters ,
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
		energy += (force_component*energy_length_parameters[i]) + (delta_length*force_length_parameters[i]) ;
		effective_spring_constant += force_length_parameters[i+1] * ( force_component + force_length_parameters[i] ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████  ██████  ██████   ██████ ███████         ██████   ██████  ██     ██ ███████ ██████  ██       █████  ██     ██ 
//  ██      ██    ██ ██   ██ ██      ██              ██   ██ ██    ██ ██     ██ ██      ██   ██ ██      ██   ██ ██     ██ 
//  █████   ██    ██ ██████  ██      █████           ██████  ██    ██ ██  █  ██ █████   ██████  ██      ███████ ██  █  ██ 
//  ██      ██    ██ ██   ██ ██      ██              ██      ██    ██ ██ ███ ██ ██      ██   ██ ██      ██   ██ ██ ███ ██ 
//  ██       ██████  ██   ██  ██████ ███████ ███████ ██       ██████   ███ ███  ███████ ██   ██ ███████ ██   ██  ███ ███  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void Spring<T,N>::spring_force_powerlaw( const T & delta_length ,
										 const std::size_t & num_force_length_parameters ,
										 const std::vector<T> & force_length_parameters ,
										 const std::vector<T> & energy_length_parameters ,
										 const std::vector<T> & stiffness_length_parameters ,
										 T & effective_spring_constant ,
										 T & force_magnitude ,
										 T & energy )
{
	effective_spring_constant = static_cast<T>(0) ;
	force_magnitude = static_cast<T>(0) ;
	energy = static_cast<T>(0) ;
	T delta_length_power ;
	for( std::size_t i = 0 ; i < num_force_length_parameters ; i+=2 ) {
		delta_length_power = std::pow( delta_length , force_length_parameters[i+1]-static_cast<T>(1) ) ;
		effective_spring_constant += stiffness_length_parameters[i] * delta_length_power ;
		delta_length_power *= delta_length ;
		force_magnitude += force_length_parameters[i] * delta_length_power ;
		energy += energy_length_parameters[i] * delta_length_power * delta_length ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████ 
//  ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██    
//  ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██    
//       ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██    
//  ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    
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
//  ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████         ██████   ██████  ██      ██    ██ ███    ██  ██████  ███    ███ ██  █████  ██      
//  ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██            ██   ██ ██    ██ ██       ██  ██  ████   ██ ██    ██ ████  ████ ██ ██   ██ ██      
//  ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██            ██████  ██    ██ ██        ████   ██ ██  ██ ██    ██ ██ ████ ██ ██ ███████ ██      
//       ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██            ██      ██    ██ ██         ██    ██  ██ ██ ██    ██ ██  ██  ██ ██ ██   ██ ██      
//  ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    ███████ ██       ██████  ███████    ██    ██   ████  ██████  ██      ██ ██ ██   ██ ███████ 
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
//  ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████         ███████ ██   ██ ██████   ██████  ███    ██ ███████ ███    ██ ████████ ██  █████  ██      
//  ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██            ██       ██ ██  ██   ██ ██    ██ ████   ██ ██      ████   ██    ██    ██ ██   ██ ██      
//  ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██            █████     ███   ██████  ██    ██ ██ ██  ██ █████   ██ ██  ██    ██    ██ ███████ ██      
//       ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██            ██       ██ ██  ██      ██    ██ ██  ██ ██ ██      ██  ██ ██    ██    ██ ██   ██ ██      
//  ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    ███████ ███████ ██   ██ ██       ██████  ██   ████ ███████ ██   ████    ██    ██ ██   ██ ███████ 
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
//  ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████         ██████  ███████ ███████ ████████         ██████   ██████  ██     ██ ███████ ██████  ██       █████  ██     ██ 
//  ██         ██    ██ ██      ██      ████   ██ ██      ██      ██              ██   ██ ██      ██         ██            ██   ██ ██    ██ ██     ██ ██      ██   ██ ██      ██   ██ ██     ██ 
//  ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████         ██████  █████   ███████    ██            ██████  ██    ██ ██  █  ██ █████   ██████  ██      ███████ ██  █  ██ 
//       ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██         ██   ██ ██           ██    ██            ██      ██    ██ ██ ███ ██ ██      ██   ██ ██      ██   ██ ██ ███ ██ 
//  ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ ███████ ██   ██ ███████ ███████    ██    ███████ ██       ██████   ███ ███  ███████ ██   ██ ███████ ██   ██  ███ ███  
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
//   ██████  ███████ ████████         ███████  ██████  █████  ██      ███████         ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████ 
//  ██       ██         ██            ██      ██      ██   ██ ██      ██              ██         ██    ██ ██      ██      ████   ██ ██      ██      ██      
//  ██   ███ █████      ██            ███████ ██      ███████ ██      █████           ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████ 
//  ██    ██ ██         ██                 ██ ██      ██   ██ ██      ██                   ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██ 
//   ██████  ███████    ██    ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ 
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
//  ██████  ███████ ███████  ██████  █████  ██      ███████         ███████ ████████ ██ ███████ ███████ ███    ██ ███████ ███████ ███████ 
//  ██   ██ ██      ██      ██      ██   ██ ██      ██              ██         ██    ██ ██      ██      ████   ██ ██      ██      ██      
//  ██████  █████   ███████ ██      ███████ ██      █████           ███████    ██    ██ █████   █████   ██ ██  ██ █████   ███████ ███████ 
//  ██   ██ ██           ██ ██      ██   ██ ██      ██                   ██    ██    ██ ██      ██      ██  ██ ██ ██           ██      ██ 
//  ██   ██ ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ███████    ██    ██ ██      ██      ██   ████ ███████ ███████ ███████ 
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

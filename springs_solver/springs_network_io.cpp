#include "springs.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef _WIN32
	#define FILESEP '\\'
#else
	#define FILESEP '/'
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   ██████ ██████  ███████  █████  ████████ ███████         ███████ ██████  ██████  ██ ███    ██  ██████          ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██          ██████  ██████       ██ 
//  ██      ██   ██ ██      ██   ██    ██    ██              ██      ██   ██ ██   ██ ██ ████   ██ ██               ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██          ██    ██ ██   ██      ██ 
//  ██      ██████  █████   ███████    ██    █████           ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████           ██    ██ ██████       ██ 
//  ██      ██   ██ ██      ██   ██    ██    ██                   ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██          ██    ██ ██   ██ ██   ██ 
//   ██████ ██   ██ ███████ ██   ██    ██    ███████ ███████ ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ ███████  ██████  ██████   █████  
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
//  ██       ██████   █████  ██████          ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██         ██████  ██ ███    ██  █████  ██████  ██    ██ 
//  ██      ██    ██ ██   ██ ██   ██         ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ████   ██ ██   ██ ██   ██  ██  ██  
//  ██      ██    ██ ███████ ██   ██         ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████           ██████  ██ ██ ██  ██ ███████ ██████    ████   
//  ██      ██    ██ ██   ██ ██   ██         ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ██  ██ ██ ██   ██ ██   ██    ██    
//  ███████  ██████  ██   ██ ██████  ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ ███████ ██████  ██ ██   ████ ██   ██ ██   ██    ██    
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
			s->force_length_parameters_tension = std::vector<T>( &read_data_T[0] , &read_data_T[NFLPT] ) ;
			delete [] read_data_T ;
			//
			std::fread( &read_data_uint32[3]  , sizeof(std::uint32_t) , 1 , file_ptr ) ;
			NFLPC = static_cast<std::size_t>(read_data_uint32[3]) ;
			read_data_T = new T [NFLPC] ;
			std::fread( read_data_T , sizeof(T) , NFLPC , file_ptr ) ;
			s->force_length_parameters_compression = std::vector<T>( &read_data_T[0] , &read_data_T[NFLPC] ) ;
			delete [] read_data_T ;
			//
			//s->precompute_parameters() ;
		}
		delete [] read_data_uint32 ;
		delete [] read_data_uint8  ;
		std::fclose( file_ptr ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ███████  █████  ██    ██ ███████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██         ██████  ██ ███    ██  █████  ██████  ██    ██ 
//  ██      ██   ██ ██    ██ ██              ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ████   ██ ██   ██ ██   ██  ██  ██  
//  ███████ ███████ ██    ██ █████           ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████           ██████  ██ ██ ██  ██ ███████ ██████    ████   
//       ██ ██   ██  ██  ██  ██              ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██          ██   ██ ██ ██  ██ ██ ██   ██ ██   ██    ██    
//  ███████ ██   ██   ████   ███████ ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ ███████ ██████  ██ ██   ████ ██   ██ ██   ██    ██    
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
			NFLPT = s->force_length_parameters_tension.size() ;
			NFLPC = s->force_length_parameters_compression.size() ;
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
// ███████  █████  ██    ██ ███████          ██████  ██    ██ ████████ ██████  ██    ██ ████████ 
// ██      ██   ██ ██    ██ ██              ██    ██ ██    ██    ██    ██   ██ ██    ██    ██    
// ███████ ███████ ██    ██ █████           ██    ██ ██    ██    ██    ██████  ██    ██    ██    
//      ██ ██   ██  ██  ██  ██              ██    ██ ██    ██    ██    ██      ██    ██    ██    
// ███████ ██   ██   ████   ███████ ███████  ██████   ██████     ██    ██       ██████     ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::save_output( void )
{
	save_network_binary( file_output_nodes.c_str() , file_output_springs.c_str() ) ;
	return ;
}

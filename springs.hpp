//
//  springs.hpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#ifndef springs_hpp
#define springs_hpp

#include <cstdlib>
#include <vector>
#include <utility>
#include <string>
#include <random>
#include <memory>

#include "vectors_nd.hpp"

///
class NetworkParameters {
public:
	//
	std::string dir_input ;
	std::string dir_output ;
	std::string file_input_parameters ;
	//
	std::size_t num_points  = 0 ;
	std::size_t num_springs = 0 ;
	std::string precision = "double" ;
	std::size_t num_dimensions = 2 ;
	std::size_t num_stiffness_tension     = 0 ;
	std::size_t num_stiffness_compression = 0 ;
	//
	NetworkParameters( const std::string & , const std::string & ) ;
	void load_parameters( const char * ) ;
	void save_parameters( const char * ) ;
} ;

///
template< class T , std::size_t N >
class Point {
public:
	Vector<T,N> position ;
	Vector<T,N> force_applied ;
	bool fixed ;
} ;

///
template< class T , std::size_t N >
class Spring {
public:
	static std::size_t num_stiffness_tension ;
	static std::size_t num_stiffness_compression ;
public:
	Point<T,N> * start ;
	Point<T,N> * end ;
	std::vector<T> stiffness_tension ;
	std::vector<T> stiffness_compression ;
	T length ;
	T rest_length ;
	Vector<T,N> force ;
	bool allow_compression ;
public:
	T spring_energy( void ) ;
	Vector<T,N> get_force( void ) { return force ; } ;
} ;

///
class ASpringNetwork {
public:
	virtual void setup( const NetworkParameters & ) = 0 ;
	virtual void solve( void ) = 0 ;
	static std::unique_ptr<ASpringNetwork> create_spring_network_obj( const NetworkParameters & ) ;
} ;

///
template< class T , std::size_t N >
class SpringNetwork : public ASpringNetwork {
	// graph data structures
	struct Link {
		Point<T,N> * point ;
		Spring<T,N> * spring ;
		T spring_direction ;
	} ;
	struct Node {
		Point<T,N> * point ;
		std::vector< Link > links ;
		Vector<T,N> net_force ;
		T net_force_magnitude = static_cast<T>(0) ;
	} ;
	// iterator types
	typedef typename std::vector<  Point<T,N> >::iterator iterPoint ;
	typedef typename std::vector< Spring<T,N> >::iterator iterSpring ;
	typedef typename std::vector<        Link >::iterator iterLink ;
	typedef typename std::vector<        Node >::iterator iterNode ;
private:
	//
	std::string dir_input ;
	std::string dir_output ;
	std::string file_input_parameters ;
	std::string file_input_nodes ;
	std::string file_input_springs ;
	std::string file_output_nodes ;
	std::string file_output_springs ;
	//
	std::size_t num_points ;
	std::size_t num_springs ;
	std::vector<  Point<T,N> > points ;
	std::vector<  Point<T,N> > points_init ;
	std::vector<  Point<T,N> > points_prev ;
	std::vector<  Point<T,N> > points_best ;
	std::vector< Spring<T,N> > springs ;
	std::vector<        Node > nodes ;
	std::default_random_engine rng ;
	std::uniform_real_distribution<T> uni_0_1 ;
	std::uniform_real_distribution<T> uni_1_1 ;
	T max_spring_length ;
public:
	//
	T total_energy( void ) ;
	void move_points_force( const T & ) ;
	void move_points_rand( const T & ) ;
	bool test_reboot( void ) ;
	bool accept_new_points( const T & , const T & ) ;
	void find_max_spring_length( void ) ;
	T heat_up( const T & , const T & ) ;
	void solve( void ) ;
	void anneal( void ) ;
	void save_output( void ) ;
	// input & output
	void setup( const NetworkParameters & ) ;
	void load_network_binary( const char * , const char * ) ;
	void save_network_binary( const char * , const char * ) ;
	void construct_network( void ) ;
} ;

#endif /* springs_hpp */

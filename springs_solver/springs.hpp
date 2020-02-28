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
#include "spmat.hpp"

///
void make_dir( const std::string & ) ;

///
class NetworkParameters {
public:
	//
	std::string dir_input ;
	std::string dir_output ;
	std::string file_input_parameters ;
	std::string file_output_parameters ;
	//
	std::size_t num_points  = 0 ;
	std::size_t num_springs = 0 ;
	std::string precision = "double" ;
	std::size_t num_dimensions = 2 ;
	std::size_t num_stiffness_tension     = 0 ;
	std::size_t num_stiffness_compression = 0 ;
	//
	std::string algorithm = "newton" ;
	std::size_t num_iter_save  = 0 ;
	std::size_t num_iter_print = 0 ;
	std::size_t num_iter_max   = 0 ;
	bool use_sum_net_force = false ;
	bool use_numerical_hessian = false ;
	double tolerance_change_energy = 1.0e-12 ;
	double tolerance_sum_net_force = 1.0e-12 ;

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
	Vector<T,N> net_force ;
	T net_force_magnitude = static_cast<T>(0) ;
	Vector<bool,N> fixed_dim ;
	bool fixed_all_dim ;
	bool not_referenced ;
public:
	void move( const Vector<T,N> & ) ;
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
	T effective_stiffness ;
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
	std::size_t num_iter_save = 0 ;
	virtual void setup( const NetworkParameters & ) = 0 ;
	virtual void solve( void ) = 0 ;
	static std::unique_ptr<ASpringNetwork> create_spring_network_obj( const NetworkParameters & ) ;
} ;

///
template< class T , std::size_t N >
class SpringNetwork : public ASpringNetwork {
	// graph data structures
	struct Link {
		std::size_t node_index ;
		Point<T,N> * point ;
		Spring<T,N> * spring ;
		T spring_direction ;
	} ;
	struct Node {
		std::size_t node_index ;
		Point<T,N> * point ;
		std::vector< Link > links ;
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
	std::string dir_output_iter ;
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
	//
	std::default_random_engine rng ;
	std::uniform_real_distribution<T> uni_0_1 ;
	std::uniform_real_distribution<T> uni_1_1 ;
	//
	T max_spring_length ;
	T scale_length    = static_cast<T>(1.0) ;
	T scale_stiffness = static_cast<T>(1.0) ;
	T scale_force     = static_cast<T>(1.0) ;
	Vector<T,N> scale_position_min ;
	Vector<T,N> scale_position_max ;
	Vector<T,N> scale_position_range ;
	//
	std::string algorithm ;
	std::size_t num_iter_save ;
	std::size_t num_iter_print ;
	std::size_t num_iter_max ;
	bool use_sum_net_force ;
	bool use_numerical_hessian ;
	T tolerance_change_energy ;
	T tolerance_sum_net_force ;
	T sum_net_force_magnitude ;
	std::vector<T> neg_gradient ;
	std::vector<T> step_direction ;
	spmat<T> hessian ;
public:
	//
	T total_energy( const bool & ) ;
	void move_points_newton( const T & ) ;
	void move_points_force( const T & ) ;
	void move_points_rand( const T & ) ;
	bool test_reboot( void ) ;
	bool accept_new_points( const T & , const T & ) ;
	void find_max_spring_length( void ) ;
	T heat_up( const T & , const T &, const bool & ) ;
	void solve( void ) ;
	void anneal( void ) ;
	void minimize_energy( void ) ;
	void minimize_energy_newton( void ) ;
	void compute_gradient( void ) ;
	void compute_hessian_numerical( void ) ;
	void compute_hessian_analytical( void ) ;
	void compute_newton_step_direction( void ) ;
	void save_output( void ) ;
	// input & output
	void setup( const NetworkParameters & ) ;
	void load_network_binary( const char * , const char * ) ;
	void save_network_binary( const char * , const char * ) ;
	void construct_network( void ) ;
	void get_scale_network( void ) ;
	void rescale_network( const int & ) ;
} ;

#endif /* springs_hpp */

//
//  springs.hpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//
// http://patorjk.com/software/taag/#p=display&h=0&v=0&c=c%2B%2B&f=ANSI%20Regular&t=function_name
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
#include "klein_summer.hpp"
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
	//
	std::string algorithm = "newton" ;
	std::string objective = "energy" ;
	std::size_t num_iter_save  = 0 ;
	std::size_t num_iter_print = 0 ;
	std::size_t num_iter_max   = 0 ;
	bool include_force_fixed_nodes = false ;
	bool use_numerical_hessian = false ;
	double tolerance_change_objective = 1.0e-12 ;
	double tolerance_sum_net_force = 1.0e-12 ;
	//
	bool parallelism_enabled = true ;
	std::size_t user_num_threads = 4 ;

	//
	NetworkParameters( const std::string & , const std::string & ) ;
	void load_parameters( const char * ) ;
	void save_parameters( const char * ) ;
	void setup_multithreading( const int ) ;
} ;

///
template< class T , std::size_t N >
class Point {
public:
	int node_index = -1 ;
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
	enum ForceLengthRelationship { none=0 , polynomial=1 , exponential=2 , powerlaw=3 } ;
public:
	Point<T,N> * start ;
	Point<T,N> * end ;
	ForceLengthRelationship force_length_type_tension ;
	ForceLengthRelationship force_length_type_compression ;
	std::size_t num_force_length_parameters_tension ;
	std::size_t num_force_length_parameters_compression ;
	std::vector<T> force_length_parameters_tension ;
	std::vector<T> force_length_parameters_compression ;
	std::vector<T> energy_length_parameters_tension ;
	std::vector<T> energy_length_parameters_compression ;
	std::vector<T> stiffness_length_parameters_tension ;
	std::vector<T> stiffness_length_parameters_compression ;
	T rest_length ;
	T effective_spring_constant ;
	T length ;
	T energy ;
	Vector<T,N> force ;
	Vector<T,N> delta_position ;
public:
	void precompute_parameters( void ) ;
	T spring_energy( void ) ;
	void spring_tension( const T & , T & , T & ) ;
	void spring_compression( const T & , T & , T & ) ;
	static void spring_force_polynomial ( const T & , const std::size_t & , const std::vector<T> & , const std::vector<T> & , const std::vector<T> & , T & , T & , T & ) ;
	static void spring_force_exponential( const T & , const std::size_t & , const std::vector<T> & , const std::vector<T> & , const std::vector<T> & , T & , T & , T & ) ;
	static void spring_force_powerlaw   ( const T & , const std::size_t & , const std::vector<T> & , const std::vector<T> & , const std::vector<T> & , T & , T & , T & ) ;
	T spring_stiffness_rest( void ) ;
	static T spring_stiffness_rest_polynomial ( const std::size_t & , const std::vector<T> & ) ;
	static T spring_stiffness_rest_exponential( const std::size_t & , const std::vector<T> & ) ;
	static T spring_stiffness_rest_powerlaw   ( const std::size_t & , const std::vector<T> & ) ;
	T spring_get_scale_stiffness( void ) ;
	void spring_rescale( const T & , const T & ) ;
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
		int node_index ;
		Point<T,N> * point ; // avoid pointer invalidation!
		Spring<T,N> * spring ; // avoid pointer invalidation!
		T spring_direction ;
	} ;
	struct Node {
		std::size_t node_index ;
		Point<T,N> * point ; // avoid pointer invalidation!
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
	std::vector<  Point<T,N>   > points ;
	std::vector<  Point<T,N>   > points_init ;
	std::vector<  Point<T,N>   > points_prev ;
	std::vector<  Point<T,N>   > points_best ;
	std::vector< Spring<T,N>   > springs ;
	std::vector<        Node   > nodes ;
	std::vector< Spring<T,N> * > springs_used ;
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
	std::string objective ;
	std::size_t num_iter_save ;
	std::size_t num_iter_print ;
	std::size_t num_iter_max ;
	bool include_force_fixed_nodes ;
	bool use_numerical_hessian ;
	T tolerance_change_objective ;
	T tolerance_sum_net_force ;
	T energy ;
	T sum_net_force_magnitude ;
	T max_net_force_magnitude ;
	std::vector<T> neg_gradient ;
	std::vector<T> step_direction ;
	spmat<T> hessian ;
	//
	KleinSummer<T> ksum ;
	//
	bool parallelism_enabled = true ;
	std::size_t num_threads ;
	std::vector< std::vector< std::pair< std::size_t , bool > > > thread_points  ; // index , shared
	std::vector< std::vector< std::pair< std::size_t , bool > > > thread_springs ; // index , responsible
	std::vector< std::vector< std::pair< Point<T,N> , std::size_t > > > points_shared ;
	std::vector< Point<T,N> * > points_shared_location ; // avoid pointer invalidation!
	//
public:
	//
	void update_springs( void ) ;
	void update_springs( std::vector< Spring<T,N>   > & ) ;
	void update_springs( std::vector< Spring<T,N> * > & ) ;
	void update_forces( void ) ;
	void update_forces( std::vector< Point<T,N>   > & , std::vector< Spring<T,N>   > & ) ;
	void update_forces( std::vector< Point<T,N> * > & , std::vector< Spring<T,N> * > & ) ;
	T total_energy( void ) ;
	T total_energy( const std::vector< Spring<T,N>   > & , const std::vector< std::pair< std::size_t , bool > > & , KleinSummer<T> & ) ;
	T total_energy( const std::vector< Spring<T,N> * > & , const std::vector< std::pair< std::size_t , bool > > & , KleinSummer<T> & ) ;
	T get_objective( void ) ;
	T get_objective( const T , const T , const T ) ;
	void get_net_force_mag( void ) ;
	void get_net_force_mag( const std::vector< Point<T,N>   > & , KleinSummer<T> & , T & , T & ) ;
	void get_net_force_mag( const std::vector< Point<T,N> * > & , KleinSummer<T> & , T & , T & ) ;
	void move_points_newton( const T & ) ;
	void move_points_force( const T & ) ;
	void move_points_force( const T & , std::vector< Point<T,N>   > & ) ;
	void move_points_force( const T & , std::vector< Point<T,N> * > & ) ;
	void move_points_rand( const T & ) ;
	void move_points_rand( const T & , std::default_random_engine & , std::uniform_real_distribution<T> & ) ;
	void move_points_rand( const T & , std::default_random_engine & , std::uniform_real_distribution<T> & , std::vector< Point<T,N>   > & ) ;
	void move_points_rand( const T & , std::default_random_engine & , std::uniform_real_distribution<T> & , std::vector< Point<T,N> * > & ) ;
	bool test_reboot( void ) ;
	bool test_reboot( const std::vector< Spring<T,N>   > & ) ;
	bool test_reboot( const std::vector< Spring<T,N> * > & ) ;
	bool accept_new_points( const T & , const T & , const T ) ;
	void find_max_spring_length( void ) ;
	T heat_up( const T & , const std::size_t & ) ;
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
	void setup_parallel( const NetworkParameters & ) ;
	void partition_domain( void ) ;
	void setup_shared_data( void ) ;
	void update_points_shared( const std::vector< std::pair< Point<T,N> * , Point<T,N> * > > & ) ;
	void update_points_all( const std::vector< Point<T,N>   > & ) ;
	void update_points_all( const std::vector< Point<T,N> * > & ) ;
	void setup_local_data( const int , std::vector<Point<T,N>> & , std::vector<Spring<T,N>> & , std::vector< std::pair< Point<T,N> * , Point<T,N> * > > & ) ;
	void load_network_binary( const char * , const char * ) ;
	void save_network_binary( const char * , const char * ) ;
	void construct_network( void ) ;
	void get_scale_network( void ) ;
	void rescale_network( const int & ) ;
} ;

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

#endif /* springs_hpp */
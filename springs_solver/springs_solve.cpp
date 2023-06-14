#include "springs.hpp"
#include "vectors_nd.hpp"
#include "spmat.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <string>
#include <random>

#include <omp.h>
#include <stdlib.h>

#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███████  ██████  ██      ██    ██ ███████ 
// ██      ██    ██ ██      ██    ██ ██      
// ███████ ██    ██ ██      ██    ██ █████   
//      ██ ██    ██ ██       ██  ██  ██      
// ███████  ██████  ███████   ████   ███████ 
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
// ████████  ██████  ████████  █████  ██              ███████ ███    ██ ███████ ██████   ██████  ██    ██ 
//    ██    ██    ██    ██    ██   ██ ██              ██      ████   ██ ██      ██   ██ ██        ██  ██  
//    ██    ██    ██    ██    ███████ ██              █████   ██ ██  ██ █████   ██████  ██   ███   ████   
//    ██    ██    ██    ██    ██   ██ ██              ██      ██  ██ ██ ██      ██   ██ ██    ██    ██    
//    ██     ██████     ██    ██   ██ ███████ ███████ ███████ ██   ████ ███████ ██   ██  ██████     ██    
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template< class T , std::size_t N >
// T SpringNetwork<T,N>::total_energy( void )
// {
// 	// using klein summation, modification of kahan summation,
// 	// to account for floating point precision issues with sums
// 	// of many small numbers, or small + large numbers
// 	ksum.reset() ;
// 	iterNode n ;
// 	iterNode n_end ;
// 	iterLink l ;
// 	iterLink l_end ;
// 	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
// 		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
// 			ksum.add( l->spring->energy ) ;
// 		}
// 	}
// 	/*
// 	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
// 		ksum.add( s->energy ) ;
// 	}
// 	//*/
// 	T energy = ksum.result() ;
// 	/*
// 	// external work done by applied loads
// 	Vector<T,N> displacement ;
// 	for( std::size_t p = 0 ; p < num_points ; ++p ) {
// 		displacement = points[p].position ;
// 		displacement -= points_init[p].position ;
// 		energy -= points[p].force_applied.dot( displacement ) ;
// 	}
// 	//*/
// 	return energy ;
// }

template< class T , std::size_t N >
T SpringNetwork<T,N>::total_energy( void )
{
	// using klein summation, modification of kahan summation,
	// to account for floating point precision issues with sums
	// of many small numbers, or small + large numbers
	ksum.reset() ;
	for( std::size_t i = 0 ; i < springs_used.size() ; ++i ) {
		ksum.add( springs_used[i]->energy ) ;
	}
	T energy = ksum.result() ;
	return energy ;
}

template< class T , std::size_t N >
T SpringNetwork<T,N>::total_energy( const std::vector< Spring<T,N> > & springs_subset ,
									const std::vector< std::pair< std::size_t , bool > > & springs_info ,
									KleinSummer<T> & ksum )
{
	// sum all spring energies in the local subset
	// some springs are shared by more than one partition, only one is responsible
	ksum.reset() ;
	for( std::size_t s = 0 ; s < springs_subset.size() ; ++s ) {
		if( springs_info[s].second ) { // responsible?
			ksum.add( springs_subset[s].energy ) ;
		}
	}
	return ksum.result() ;
}

template< class T , std::size_t N >
T SpringNetwork<T,N>::total_energy( const std::vector< Spring<T,N> * > & springs_subset ,
									const std::vector< std::pair< std::size_t , bool > > & springs_info ,
									KleinSummer<T> & ksum )
{
	// sum all spring energies in the local subset
	// some springs are shared by more than one partition, only one is responsible
	ksum.reset() ;
	for( std::size_t s = 0 ; s < springs_subset.size() ; ++s ) {
		if( springs_info[s].second ) { // responsible?
			ksum.add( springs_subset[s]->energy ) ;
		}
	}
	return ksum.result() ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██    ██ ██████  ██████   █████  ████████ ███████         ███████ ██████  ██████  ██ ███    ██  ██████  ███████ 
// ██    ██ ██   ██ ██   ██ ██   ██    ██    ██              ██      ██   ██ ██   ██ ██ ████   ██ ██       ██      
// ██    ██ ██████  ██   ██ ███████    ██    █████           ███████ ██████  ██████  ██ ██ ██  ██ ██   ███ ███████ 
// ██    ██ ██      ██   ██ ██   ██    ██    ██                   ██ ██      ██   ██ ██ ██  ██ ██ ██    ██      ██ 
//  ██████  ██      ██████  ██   ██    ██    ███████ ███████ ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::update_springs( void )
{
	// internal energy due to spring tension/compression
	// parallelize computation for more than 10^4 springs
	for( std::size_t i = 0 ; i < springs_used.size() ; ++i ) {
		springs_used[i]->spring_energy() ;
	}
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::update_springs( std::vector< Spring<T,N> > & springs_subset )
{
	// internal energy due to spring tension/compression
	for( std::size_t i = 0 ; i < springs_subset.size() ; ++i ) {
		springs_subset[i].spring_energy() ;
	}
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::update_springs( std::vector< Spring<T,N> * > & springs_subset )
{
	// internal energy due to spring tension/compression
	for( std::size_t i = 0 ; i < springs_subset.size() ; ++i ) {
		springs_subset[i]->spring_energy() ;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██    ██ ██████  ██████   █████  ████████ ███████         ███████  ██████  ██████   ██████ ███████ ███████ 
// ██    ██ ██   ██ ██   ██ ██   ██    ██    ██              ██      ██    ██ ██   ██ ██      ██      ██      
// ██    ██ ██████  ██   ██ ███████    ██    █████           █████   ██    ██ ██████  ██      █████   ███████ 
// ██    ██ ██      ██   ██ ██   ██    ██    ██              ██      ██    ██ ██   ██ ██      ██           ██ 
//  ██████  ██      ██████  ██   ██    ██    ███████ ███████ ██       ██████  ██   ██  ██████ ███████ ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template< class T , std::size_t N >
// void SpringNetwork<T,N>::update_forces( void )
// {
// 	/*
// 	// account for force imbalances, "potential energy"
// 	// combine external applied force & internal spring forces
// 	iterNode n ;
// 	iterNode n_end ;
// 	iterLink l ;
// 	iterLink l_end ;
// 	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
// 		n->point->net_force = n->point->force_applied ;
// 		for( l = n->links.begin() , l_end = n->links.end() ; l != l_end ; ++l ) {
// 			if( l->spring_direction > static_cast<T>(0) ) {
// 				n->point->net_force += l->spring->get_force() ;
// 			} else {
// 				n->point->net_force -= l->spring->get_force() ;
// 			}
// 		}
// 		n->point->net_force_magnitude = n->point->net_force.norm() ;
// 	}
// 	sum_net_force_magnitude = static_cast<T>(0) ;
// 	for( iterNode n = nodes.begin() ; n != nodes.end() ; ++n ) {
// 		sum_net_force_magnitude += n->point->net_force_magnitude ;
// 	}
// 	//*/

// 	// net force on each point (including fixed points!)
// 	/*
// 	for( std::size_t p = 0 ; p < num_points ; ++p ) {
// 		points[p].net_force = static_cast<T>(0.0) ;
// 	}
// 	//*/
// 	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
// 		p->net_force = static_cast<T>(0.0) ;
// 	}
// 	for( std::size_t i = 0 ; i < springs_used.size() ; ++i ) {
// 		springs_used[i]->start->net_force += springs_used[i]->force ;
// 		springs_used[i]->end  ->net_force -= springs_used[i]->force ;
// 	}
	
// 	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
// 		s->start->net_force += s->force ;
// 		s->end->net_force   -= s->force ;
// 	}
// 	//
// 	for( std::size_t p = 0 ; p < num_points ; ++p ) {
// 		Point<T,N> * current_point = &points[p] ;
// 		current_point->net_force += current_point->force_applied ;
// 		current_point->net_force_magnitude = current_point->net_force.norm() ;
// 	}
// 	/*
// 	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
// 		p->net_force = static_cast<T>(0.0) ;
// 	}
// 	for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
// 		Vector<T,N> spring_force = s->get_force() ;
// 		s->start->net_force += spring_force ;
// 		s->end->net_force   -= spring_force ;
// 	}
// 	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
// 		p->net_force += p->force_applied ;
// 		p->net_force_magnitude = p->net_force.norm() ;
// 	}
// 	//*/
// 	ksum.reset() ;
// 	max_net_force_magnitude = static_cast<T>(0.0) ;
// 	if( include_force_fixed_nodes ) {
// 		for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
// 			if( !p->not_referenced ) {
// 				max_net_force_magnitude = ( p->net_force_magnitude > max_net_force_magnitude ) ? p->net_force_magnitude : max_net_force_magnitude ;
// 				ksum.add( p->net_force_magnitude ) ;
// 			}
// 		}
// 	} else {
// 		for( iterNode n = nodes.begin() ; n != nodes.end() ; ++n ) {
// 			max_net_force_magnitude = ( n->point->net_force_magnitude > max_net_force_magnitude ) ? n->point->net_force_magnitude : max_net_force_magnitude ;
// 			ksum.add( n->point->net_force_magnitude ) ;
// 		}
// 	}
// 	sum_net_force_magnitude = ksum.result() ;
// 	return ;
// }

template< class T , std::size_t N >
void SpringNetwork<T,N>::update_forces( void )
{
	// net force on each point (including fixed points!)
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		p->net_force = static_cast<T>(0.0) ;
	}
	for( std::size_t i = 0 ; i < springs_used.size() ; ++i ) {
		springs_used[i]->start->net_force += springs_used[i]->force ;
		springs_used[i]->end  ->net_force -= springs_used[i]->force ;
	}
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		p->net_force += p->force_applied ;
		p->net_force_magnitude = p->net_force.norm() ;
	}
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::update_forces( std::vector< Point<T,N> > & points_subset ,
										std::vector< Spring<T,N> > & springs_subset )
{
	// net force on each point (including fixed points!)
	for( iterPoint p = points_subset.begin() ; p != points_subset.end() ; ++p ) {
		p->net_force = static_cast<T>(0.0) ;
	}
	for( iterSpring s = springs_subset.begin() ; s != springs_subset.end() ; ++s ) {
		s->start->net_force += s->force ;
		s->end->net_force   -= s->force ;
	}
	for( iterPoint p = points_subset.begin() ; p != points_subset.end() ; ++p ) {
		p->net_force += p->force_applied ;
		p->net_force_magnitude = p->net_force.norm() ;
	}
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::update_forces( std::vector< Point<T,N> * > & points_subset ,
										std::vector< Spring<T,N> * > & springs_subset )
{
	// net force on each point (including fixed points!)
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		points_subset[p]->net_force = static_cast<T>(0.0) ;
	}
	for( std::size_t s = 0 ; s < springs_subset.size() ; ++s ) {
		springs_subset[s]->start->net_force += springs_subset[s]->force ;
		springs_subset[s]->end->net_force   -= springs_subset[s]->force ;
	}
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		points_subset[p]->net_force += points_subset[p]->force_applied ;
		points_subset[p]->net_force_magnitude = points_subset[p]->net_force.norm() ;
	}
	return ;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   ██████  ███████ ████████         ███    ██ ███████ ████████         ███████  ██████  ██████   ██████ ███████         ███    ███  █████   ██████  
//  ██       ██         ██            ████   ██ ██         ██            ██      ██    ██ ██   ██ ██      ██              ████  ████ ██   ██ ██       
//  ██   ███ █████      ██            ██ ██  ██ █████      ██            █████   ██    ██ ██████  ██      █████           ██ ████ ██ ███████ ██   ███ 
//  ██    ██ ██         ██            ██  ██ ██ ██         ██            ██      ██    ██ ██   ██ ██      ██              ██  ██  ██ ██   ██ ██    ██ 
//   ██████  ███████    ██    ███████ ██   ████ ███████    ██    ███████ ██       ██████  ██   ██  ██████ ███████ ███████ ██      ██ ██   ██  ██████  
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::get_net_force_mag( void )
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
		//*
		for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
			if( !p->not_referenced && !p->fixed_any_dim ) {
				max_net_force_magnitude = ( p->net_force_magnitude > max_net_force_magnitude ) ? p->net_force_magnitude : max_net_force_magnitude ;
				ksum.add( p->net_force_magnitude ) ;
			}
		}
		//*/
		/*
		for( iterNode n = nodes.begin() ; n != nodes.end() ; ++n ) {
			max_net_force_magnitude = ( n->point->net_force_magnitude > max_net_force_magnitude ) ? n->point->net_force_magnitude : max_net_force_magnitude ;
			ksum.add( n->point->net_force_magnitude ) ;
		}
		//*/
	}
	sum_net_force_magnitude = ksum.result() ;
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::get_net_force_mag( const std::vector< Point<T,N> > & points_subset ,
											KleinSummer<T> & ksum ,
											T & max_net_force_magnitude ,
											T & sum_net_force_magnitude )
{
	ksum.reset() ;
	max_net_force_magnitude = static_cast<T>(0.0) ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		max_net_force_magnitude = ( points_subset[p].net_force_magnitude > max_net_force_magnitude ) ? points_subset[p].net_force_magnitude : max_net_force_magnitude ;
		ksum.add( points_subset[p].net_force_magnitude ) ;
	}
	sum_net_force_magnitude = ksum.result() ;
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::get_net_force_mag( const std::vector< Point<T,N> * > & points_subset ,
											KleinSummer<T> & ksum ,
											T & max_net_force_magnitude ,
											T & sum_net_force_magnitude )
{
	ksum.reset() ;
	max_net_force_magnitude = static_cast<T>(0.0) ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		max_net_force_magnitude = ( points_subset[p]->net_force_magnitude > max_net_force_magnitude ) ? points_subset[p]->net_force_magnitude : max_net_force_magnitude ;
		ksum.add( points_subset[p]->net_force_magnitude ) ;
	}
	sum_net_force_magnitude = ksum.result() ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ███████ ████████          ██████  ██████       ██ ███████  ██████ ████████ ██ ██    ██ ███████ 
// ██       ██         ██            ██    ██ ██   ██      ██ ██      ██         ██    ██ ██    ██ ██      
// ██   ███ █████      ██            ██    ██ ██████       ██ █████   ██         ██    ██ ██    ██ █████   
// ██    ██ ██         ██            ██    ██ ██   ██ ██   ██ ██      ██         ██    ██  ██  ██  ██      
//  ██████  ███████    ██    ███████  ██████  ██████   █████  ███████  ██████    ██    ██   ████   ███████ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
T SpringNetwork<T,N>::get_objective( void )
{
	if( objective=="energy"   ) { return energy                  ; } else
	if( objective=="sumforce" ) { return sum_net_force_magnitude ; } else
	if( objective=="maxforce" ) { return max_net_force_magnitude ; }
	else { return static_cast<T>(0.0) ; }
}

template< class T , std::size_t N >
T SpringNetwork<T,N>::get_objective( const T energy , const T sum_net_force_magnitude , const T max_net_force_magnitude )
{
	if( objective=="energy"   ) { return energy                  ; } else
	if( objective=="sumforce" ) { return sum_net_force_magnitude ; } else
	if( objective=="maxforce" ) { return max_net_force_magnitude ; }
	else { return static_cast<T>(0.0) ; }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███    ███  ██████  ██    ██ ███████         ██████   ██████  ██ ███    ██ ████████ ███████         ███████  ██████  ██████   ██████ ███████ 
// ████  ████ ██    ██ ██    ██ ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ██      ██    ██ ██   ██ ██      ██      
// ██ ████ ██ ██    ██ ██    ██ █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         █████   ██    ██ ██████  ██      █████   
// ██  ██  ██ ██    ██  ██  ██  ██              ██      ██    ██ ██ ██  ██ ██    ██         ██         ██      ██    ██ ██   ██ ██      ██      
// ██      ██  ██████    ████   ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ██       ██████  ██   ██  ██████ ███████ 
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
	//#pragma omp for schedule(static)
	for( std::size_t n = 0 ; n < nodes.size() ; ++n ) {
		Vector<T,N> small_displacement ;
		small_displacement = nodes[n].point->net_force ;
		small_displacement *= step_size ;
		nodes[n].point->move( small_displacement ) ;
	}
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_force( const T & step_size , std::vector< Point<T,N> > & points_subset )
{
	// apply net force & move small displacement towards equilibrating position
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		Vector<T,N> small_displacement ;
		small_displacement = points_subset[p].net_force ;
		small_displacement *= step_size ;
		points_subset[p].move( small_displacement ) ;
	}
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_force( const T & step_size , std::vector< Point<T,N> * > & points_subset )
{
	// apply net force & move small displacement towards equilibrating position
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		Vector<T,N> small_displacement ;
		small_displacement = points_subset[p]->net_force ;
		small_displacement *= step_size ;
		points_subset[p]->move( small_displacement ) ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███    ███  ██████  ██    ██ ███████         ██████   ██████  ██ ███    ██ ████████ ███████         ██████   █████  ███    ██ ██████  
// ████  ████ ██    ██ ██    ██ ██              ██   ██ ██    ██ ██ ████   ██    ██    ██              ██   ██ ██   ██ ████   ██ ██   ██ 
// ██ ████ ██ ██    ██ ██    ██ █████           ██████  ██    ██ ██ ██ ██  ██    ██    ███████         ██████  ███████ ██ ██  ██ ██   ██ 
// ██  ██  ██ ██    ██  ██  ██  ██              ██      ██    ██ ██ ██  ██ ██    ██         ██         ██   ██ ██   ██ ██  ██ ██ ██   ██ 
// ██      ██  ██████    ████   ███████ ███████ ██       ██████  ██ ██   ████    ██    ███████ ███████ ██   ██ ██   ██ ██   ████ ██████  
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

template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_rand( const T & amplitude ,
										   std::default_random_engine & rng ,
										   std::uniform_real_distribution<T> & uni )
{
	iterNode n ;
	iterNode n_end ;
	Vector<T,N> small_displacement ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( std::size_t d = 0 ; d < N ; ++d ) {
			small_displacement[d] = amplitude * uni(rng) ;
		}
		n->point->move( small_displacement ) ;
	}
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_rand( const T & amplitude ,
										   std::default_random_engine & rng ,
										   std::uniform_real_distribution<T> & uni ,
										   std::vector< Point<T,N> > & points_subset )
{
	Vector<T,N> small_displacement ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		for( std::size_t d = 0 ; d < N ; ++d ) {
			small_displacement[d] = amplitude * uni(rng) ;
		}
		points_subset[p].move( small_displacement ) ;
	}
	return ;
}

template< class T , std::size_t N >
void SpringNetwork<T,N>::move_points_rand( const T & amplitude ,
										   std::default_random_engine & rng ,
										   std::uniform_real_distribution<T> & uni ,
										   std::vector< Point<T,N> * > & points_subset )
{
	Vector<T,N> small_displacement ;
	for( std::size_t p = 0 ; p < points_subset.size() ; ++p ) {
		for( std::size_t d = 0 ; d < N ; ++d ) {
			small_displacement[d] = amplitude * uni(rng) ;
		}
		points_subset[p]->move( small_displacement ) ;
	}
	return ;
}





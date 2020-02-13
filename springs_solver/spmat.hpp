//
//  spmat.hpp
//  HW4
//
//  Created by Jacob Herrmann on 4/5/16.
//  Copyright © 2016 Jacob Herrmann. All rights reserved.
//

#ifndef spmat_h
#define spmat_h

#include <vector>
#include <utility>
#include <map>

template< class T >
class spmat {
	
	// nested map for sparse matrix container
	// using Compressed-Row Storage:
	//   * outer map indexes rows
	//   * inner map indexes columns
	// this strategy is ideal for:
	//   * quickly indexing large arrays
	//   * performing insertions and deletions (e.g. boundary conditions)
	typedef          std::vector< std::vector< std::pair< std::size_t , T > > > mat_T ;
	typedef typename std::vector< std::vector< std::pair< std::size_t , T > > >::iterator itRow ;
	typedef typename std::vector< std::pair< std::size_t , T > >::iterator itCol ;
	typedef          std::vector< T > vec_T ;
	typedef typename std::vector< T >::iterator itVec ;
	
public:
	
	 /**************************************************/
	/*********** DEFINED IN SPMAT_CLASS.CPP ***********/

	// constructor
	spmat() ;
	spmat( std::size_t numRow ,   // number of rows
		   std::size_t numCol ) ; // number of columns
	
	// destructor
	~spmat( void ) ;
	
	// access matrix values
	T get_val( std::size_t indRow ,   // row index (must be less than numRow)
			   std::size_t indCol ) ; // column index (must be less than numCol)
	
	// modify matrix values
	void set_val( std::size_t indRow , // row index (must be less than numRow)
				  std::size_t indCol , // column index (must be less than numCol
				  T val ) ;            // new value for specificied element
	void set_val( std::vector<std::size_t> indRow , // row indexes (must be less than numRow)
				  std::vector<std::size_t> indCol , // column indexes (must be less than numCol
				  vec_T val ) ;                     // new values for specificied element
	void set_val( std::size_t *indRow ,  // row indexes (must be less than numRow)
				  std::size_t *indCol ,  // column indexes (must be less than numCol
				  T *val ,               // new values for specificied element
				  std::size_t numVal ) ; // number of elements in arrays
	
	// print matrix values
	void print_val( void ) ;
	
	// matrix-vector multiplication
	vec_T operator*( const vec_T &x ) ; // vector to be multiplied by matrix
	
	 /**************************************************/
	/********** DEFINED IN SPMAT_SOLVERS.CPP **********/
	
	// create modified matrix for more efficient solves
	void prepare_SOR() ;
	
	// solver: Successive Over-Relaxation
	void solve_SOR(
		const vec_T & b ,           // source vector
		vec_T & x0 ,                // initial guess, and solution target
		const T errTol ,            // error tolerance, convergence criterion
		const std::size_t maxIter , // maximum number of iterations
		const bool modifyOmega ) ;  // dynamically modify relaxation parameter
	
	 /**************************************************/
	/**************************************************/
	
	// accessors
	std::size_t get_numRow( void ) { return numRow ; }
	std::size_t get_numCol( void ) { return numCol ; }
	std::size_t get_numNZ ( void ) { return numNZ  ; }
	
private:
	std::size_t numRow    ; // number of rows
	std::size_t numCol    ; // number of columns
	std::size_t numNZ     ; // number of nonzeros
	mat_T M               ; // matrix container
	vec_T M_diag          ; // diagonal, for modifying matrix
	mat_T M_mod           ; // modified matrix
	std::vector<vec_T> LU ; // LU-decomposed matrix container (dense)
	
} ;

#endif // spmat_h

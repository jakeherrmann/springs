//
//  spmat_solvers.cpp
//  HW4
//
//  Created by Jacob Herrmann on 4/5/16.
//  Copyright © 2016 Jacob Herrmann. All rights reserved.
//

#include "spmat.hpp"
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <iostream>

/** create modified matrix for more efficient solves **/
template< class T >
void spmat<T>::prepare_SOR()
{
	// MODIFY COEFFICIENT MATRIX FOR MORE EFFICIENT SOLVE
	M_mod = M ;
	M_diag = vec_T( numRow , static_cast<T>(0) ) ;
	
	// (*cc).first  contains the column index
	// (*cc).second contains the nonzero element value
	itCol cc ;
	std::size_t jj ;
	
	for( jj = 0 ; jj < M_mod.size() ; ++jj ) {
		for( cc = M_mod[jj].begin() ; cc != M_mod[jj].end() ; ++cc ) {
			if( (*cc).first == jj ) {
				M_diag[jj] = (*cc).second ;
				M_mod[jj].erase( cc ) ;
				break ;
			}
		}
		for( cc = M_mod[jj].begin() ; cc != M_mod[jj].end() ; ++cc ) {
			(*cc).second /= M_diag[jj] ;
		}
	}
	
	return ;
}

/** solver: Successive Over-Relaxation **/
template< class T >
void spmat<T>::solve_SOR(
	const vec_T & b ,           // source vector
	vec_T & x0 ,                // initial guess, and solution target
	const T errTol ,            // error tolerance, convergence criterion
	const std::size_t maxIter , // maximum number of iterations
	const bool modifyOmega )    // dynamically modify relaxation parameter
{
	// output vector, initiallized with size (numRow,1) and all zeros
	vec_T x1( numRow , static_cast<T>(0) ) ;
	
	// (*rr)  contains the rows
	itRow rr ;
	
	// (*cc).first  contains the column index
	// (*cc).second contains the nonzero element value
	itCol cc ;
	
	// iterator for vector<T>
	itVec bb ;
	itVec xx ;
	
	std::size_t ii ;
	std::size_t jj ;
	
	// relaxation parameter
	T omega    = static_cast<T>(0.2) ;
	T omegaMin = static_cast<T>(0.2) ;
	T omegaMax = static_cast<T>(1.8) ;
	
	// residual, error
	T err = static_cast<T>(0) ;
	T errPrev = static_cast<T>(0) ;
	
	// temporary variables for loops
	T tmp_T ;
	
	// MODFIY COEFFICIENT MATRIX FOR MORE EFFICIENT SOLVE
	if( M_diag.empty() ) {
		prepare_SOR() ;
	}
	vec_T b_mod( numRow , static_cast<T>(0) ) ;
	for( jj = 0 ; jj < numRow ; ++jj ) {
		b_mod[jj] = b[jj] / M_diag[jj] ;
	}

	//
	for( ii = 0 ; ii < maxIter ; ++ii ) {
		
		// THIS PART IS SLOOOOOOOOW
		// find new guess using Gauss-Seidel iteration
		x1 = x0 ;
		for( rr = M_mod.begin() , bb = b_mod.begin() , xx = x1.begin() ; rr != M_mod.end() ; ++rr , ++bb , ++xx ) {
			tmp_T = (*bb) ;
			for( cc = (*rr).begin() ; cc != (*rr).end() ; ++cc ) {
				tmp_T -= (*cc).second * x1[ (*cc).first ] ;
			}
			(*xx) = tmp_T ;
		}
		
		// modify new guess using relaxation parameter
		// compute change in guess (err)
		err = static_cast<T>(0) ;
		for( jj = 0 ; jj < numRow ; ++jj ) {
			x1[jj] = x0[jj] + ( omega * ( x1[jj] - x0[jj] ) ) ;
			tmp_T = std::abs( x0[jj] - x1[jj] ) ;
			//err = ( tmp_T > err ) ? tmp_T : err ; // maximum error
			err += tmp_T ; // sum error
		}
		err /= static_cast<T>( numRow ) ;
		
		// update guess and check convergence
		x0 = x1 ;
		if( err < errTol ) { break ; }
		
		// dynamically modify relaxation parameter
		if( modifyOmega && (ii>0) ) {
			omega = ( errPrev / err ) ;
			omega = std::max( omegaMin , std::min( omegaMax , omega ) ) ;
			errPrev = err ;
		}
		
		///std::printf( "\n\t\t\t%u iter - %E err - %0.1f omega" , static_cast<int>(ii) , err , omega ) ;
	}
	
	///std::printf( "\n\t\t%u iterations - %E error" , static_cast<int>(ii) , err ) ;
	
	return ;
}


/** EXPLICIT TEMPLATE INSTANTIATIONS **/
template class spmat<float > ;
template class spmat<double> ;


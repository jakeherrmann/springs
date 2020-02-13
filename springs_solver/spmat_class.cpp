//
//  spmat.cpp
//  HW4
//
//  Created by Jacob Herrmann on 2/9/16.
//  Copyright © 2016 Jacob Herrmann. All rights reserved.
//

#include "spmat.hpp"
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <utility>
#include <algorithm>

/** constructor **/
template< class T >
spmat<T>::spmat()
{
	return ;
}

template< class T >
spmat<T>::spmat( std::size_t numRow , // number of rows
				 std::size_t numCol ) // number of columns
{
	this->numRow = numRow ;
	this->numCol = numCol ;
	this->M = typename spmat<T>::mat_T( numRow ) ;
	return ;
}

/** destructor **/
template< class T >
spmat<T>::~spmat( void )
{
	return ;
}

/** accessor for values **/
template< class T >
T spmat<T>::get_val( std::size_t indRow , // row index (must be less than numRow)
					 std::size_t indCol ) // column index (must be less than numCol)
{
	spmat<T>::itCol cc ;
	for( cc = M[indRow].begin() ; cc != M[indRow].end() ; ++cc ) {
		if( (*cc).first == indCol ) {
			return (*cc).second ;
		}
	}
	return static_cast<T>(0) ;
}

/** mutator for values **/
template< class T >
void spmat<T>::set_val( std::size_t indRow , // row index (must be less than numRow)
					    std::size_t indCol , // column index (must be less than numCol)
					    T val )              // new value for specificied element
{
	spmat<T>::itCol cc ;
	if( val != static_cast<T>(0) ) {
		for( cc = M[indRow].begin() ; cc != M[indRow].end() ; ++cc ) {
			if( (*cc).first == indCol ) {
				// adjust existing nonzero value
				(*cc).second = val ;
				break ;
			}
		}
		if( cc == M[indRow].end() ) {
			// insert nonzero value into specified element
			M[indRow].push_back( std::make_pair( indCol , val ) ) ;
		}
	} else {
		for( cc = M[indRow].begin() ; cc != M[indRow].end() ; ++cc ) {
			if( (*cc).first == indCol ) {
				// remove zero element, if it was previously nonzero
				M[indRow].erase( cc ) ;
				break ;
			}
		}
	}
	return ;
}
template< class T >
void spmat<T>::set_val( std::vector<std::size_t> indRow , // row indexes (must be less than numRow)
					    std::vector<std::size_t> indCol , // column indexes (must be less than numCol)
					    vec_T val )                       // new values for specificied element
{
	spmat<T>::itCol cc ;
	for( int jj = 0 ; jj < val.size() ; ++jj ) {
		if( val[jj] != static_cast<T>(0) ) {
			for( cc = M[ indRow[jj] ].begin() ; cc != M[ indRow[jj] ].end() ; ++cc ) {
				if( (*cc).first == indCol[jj] ) {
					// adjust existing nonzero value
					(*cc).second = val[jj] ;
					break ;
				}
			}
			if( cc == M[ indRow[jj] ].end() ) {
				// insert nonzero value into specified element
				M[ indRow[jj] ].push_back( std::make_pair( indCol[jj] , val[jj] ) ) ;
			}
		} else {
			for( cc = M[ indRow[jj] ].begin() ; cc != M[ indRow[jj] ].end() ; ++cc ) {
				if( (*cc).first == indCol[jj] ) {
					// remove zero element, if it was previously nonzero
					M[ indRow[jj] ].erase( cc ) ;
					break ;
				}
			}
		}
	}
	return ;
}
template< class T >
void spmat<T>::set_val( std::size_t *indRow , // row indexes (must be less than numRow)
					    std::size_t *indCol , // column indexes (must be less than numCol)
					    T *val ,              // new values for specificied element
					    std::size_t numVal )  // number of elements in arrays
{
	spmat<T>::itCol cc ;
	for( int jj = 0 ; jj < numVal ; ++jj ) {
		if( val[jj] != static_cast<T>(0) ) {
			for( cc = M[ indRow[jj] ].begin() ; cc != M[ indRow[jj] ].end() ; ++cc ) {
				if( (*cc).first == indCol[jj] ) {
					// adjust existing nonzero value
					(*cc).second = val[jj] ;
					break ;
				}
			}
			if( cc == M[ indRow[jj] ].end() ) {
				// insert nonzero value into specified element
				M[ indRow[jj] ].push_back( std::make_pair( indCol[jj] , val[jj] ) ) ;
			}
		} else {
			for( cc = M[ indRow[jj] ].begin() ; cc != M[ indRow[jj] ].end() ; ++cc ) {
				if( (*cc).first == indCol[jj] ) {
					// remove zero element, if it was previously nonzero
					M[ indRow[jj] ].erase( cc ) ;
					break ;
				}
			}
		}
	}
	return ;
}

/** print matrix values **/
template< class T >
void spmat<T>::print_val( void )
{
	T val ;
	for( std::size_t rr = 0 ; rr < numRow ; ++rr ) {
		for( std::size_t cc = 0 ; cc < numCol ; ++cc ) {
			val = spmat<T>::get_val( rr , cc ) ;
			if( val != 0 ) {
				std::printf( "\t%+7.2f" , val ) ;
			} else {
				std::printf( "\t   0.  " ) ;
			}
		}
		std::printf( "\n" ) ;
	}
	return ;
}

/** matrix-vector multiplication **/
template< class T >
typename spmat<T>::vec_T spmat<T>::operator*( const vec_T &x ) // vector to be multiplied by matrix
{
	// output vector, initiallized with size (numRow,1) and all zeros
	vec_T y( numRow , static_cast<T>(0) ) ;
	
	// (*rr) contains the first row
	itRow rr ;
	itVec vv ;
	
	// (*cc).first  contains the column index
	// (*cc).second contains the nonzero element value
	itCol cc ;
	
	// iterate through nonzero elements and perform each dot product
	for( rr = M.begin() , vv = y.begin() ; rr != M.end() ; ++rr , ++vv ) {
		for( cc = (*rr).begin() ; cc != (*rr).end() ; ++cc ) {
			(*vv) += (*cc).second * x[ (*cc).first ] ;
		}
	}
	return y ;
}

/** EXPLICIT TEMPLATE INSTANTIATIONS **/
template class spmat<float > ;
template class spmat<double> ;


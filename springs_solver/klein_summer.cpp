//
//  klein_summer.cpp
//  springs
//
//  Created by Jacob Herrmann on 22/03/2020.
//  Copyright © 2020 Jacob Herrmann. All rights reserved.
//

#include "klein_summer.hpp"

#include <cmath>

template< class T >
KleinSummer<T>::KleinSummer()
{
	s   = static_cast<T>(0.0) ;
	cs  = static_cast<T>(0.0) ;
	ccs = static_cast<T>(0.0) ;
	return ;
}
template< class T >
KleinSummer<T>::KleinSummer( const T init )
{
	s   = init ;
	cs  = static_cast<T>(0.0) ;
	ccs = static_cast<T>(0.0) ;
	return ;
}

template< class T >
void KleinSummer<T>::reset( void )
{
	s   = static_cast<T>(0.0) ;
	cs  = static_cast<T>(0.0) ;
	ccs = static_cast<T>(0.0) ;
	return ;
}
template< class T >
void KleinSummer<T>::reset( const T init )
{
	s   = init ;
	cs  = static_cast<T>(0.0) ;
	ccs = static_cast<T>(0.0) ;
	return ;
}

template< class T >
void KleinSummer<T>::add( const T next )
{
	T t , c , ct , cc ;
	t = s + next ;
	if( std::fabs(s) >= std::fabs(next) ) {
		c = (s-t) + next ;
	} else {
		c = (next-t) + s ;
	}
	s = t ;
	ct = cs + c ;
	if( std::fabs(cs) >= std::fabs(c) ) {
		cc = (cs-ct) + c ;
	} else {
		cc = (c-ct) + cs ;
	}
	cs = ct ;
	ccs = ccs + cc ;
	return ;
}

template< class T >
T KleinSummer<T>::result( void )
{
	return s + cs + ccs ;
}

/// EXPLICIT TEMPLATE INSTANTIATIONS ///
template class KleinSummer<float > ;
template class KleinSummer<double> ;
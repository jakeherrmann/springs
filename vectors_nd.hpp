//
//  vectors_nd.hpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#ifndef vectors_nd_hpp
#define vectors_nd_hpp

#include <array>

template< class T , std::size_t N >
class Vector {

private:
	std::array<T,N> x ;
	
public:
	// access elements
	T   operator[]( const std::size_t & ) const ; // const, cannot change indexed object
	T & operator[]( const std::size_t & ) ; // non-const, can change indexed object
	
	// scalar assignment
	Vector<T,N> & operator= ( const T & ) ;
	
	//
	T dot( const Vector<T,N> & ) const ;
	T distance( const Vector<T,N> & ) const ;
	T norm( void ) const ;
	
	// element-wise vector operations
	Vector<T,N>   operator+ ( const Vector<T,N> & ) const ;
	Vector<T,N> & operator+=( const Vector<T,N> & ) ;
	Vector<T,N>   operator- ( const Vector<T,N> & ) const ;
	Vector<T,N> & operator-=( const Vector<T,N> & ) ;
	Vector<T,N>   operator* ( const Vector<T,N> & ) const ;
	Vector<T,N> & operator*=( const Vector<T,N> & ) ;
	Vector<T,N>   operator/ ( const Vector<T,N> & ) const ;
	Vector<T,N> & operator/=( const Vector<T,N> & ) ;
	
	// scalar operations
	Vector<T,N>   operator+ ( const T & ) const ;
	Vector<T,N> & operator+=( const T & ) ;
	Vector<T,N>   operator- ( const T & ) const ;
	Vector<T,N> & operator-=( const T & ) ;
	Vector<T,N>   operator* ( const T & ) const ;
	Vector<T,N> & operator*=( const T & ) ;
	Vector<T,N>   operator/ ( const T & ) const ;
	Vector<T,N> & operator/=( const T & ) ;
} ;

#endif /* vectors_nd_hpp */

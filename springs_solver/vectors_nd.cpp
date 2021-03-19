//
//  vectors_nd.cpp
//  springs
//
//  Created by Jacob Herrmann on 13/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#include "vectors_nd.hpp"

#include <cmath>

// indexing
template< class T , std::size_t N >
inline T Vector<T,N>::operator[]( const std::size_t index )
const
{
	return x[index] ;
}

template< class T , std::size_t N >
inline T & Vector<T,N>::operator[]( const std::size_t index )
{
	return x[index] ;
}

// scalar assignment
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator=( const T a )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] = a ;
	}
	return *this ;
}

// elementwise boolean multiplication
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::zero_where( const Vector<bool,N> & y )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = ( y[n] ) ? static_cast<T>(0.0) : x[n] ;
	}
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::zero_where_not( const Vector<bool,N> & y )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = ( y[n] ) ? x[n] : static_cast<T>(0.0) ;
	}
	return v ;
}

// dot product of two vectors
template< class T , std::size_t N >
inline T Vector<T,N>::dot( const Vector<T,N> & y )
const
{
	T s = static_cast<T>(0) ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		s += x[n] * y[n] ;
	}
	return s ;
}

// distance between two points
template< class T , std::size_t N >
inline T Vector<T,N>::distance( const Vector<T,N> & y )
const
{
	T d ;
	T s = static_cast<T>(0) ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		d = x[n] - y[n] ;
		s += d*d ;
	}
	return static_cast<T>(std::sqrt(s)) ;
}

// vector norm
template< class T , std::size_t N >
inline T Vector<T,N>::norm( void )
const
{
	T s = static_cast<T>(0) ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		s += x[n]*x[n] ;
	}
	return static_cast<T>(std::sqrt(s)) ;
}

// max and min
template< class T , std::size_t N >
inline T Vector<T,N>::max( void )
const
{
	T s = x[0] ;
	for( std::size_t n = 1 ; n < N ; ++n ) {
		s = ( x[n] > s ) ? x[n] : s ;
	}
	return s ;
}
template< class T , std::size_t N >
inline T Vector<T,N>::min( void )
const
{
	T s = x[0] ;
	for( std::size_t n = 1 ; n < N ; ++n ) {
		s = ( x[n] < s ) ? x[n] : s ;
	}
	return s ;
}


// addition of two vectors
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator+( const Vector<T,N> & y )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = x[n] + y[n] ;
	}
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator+=( const Vector<T,N> & y )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] += y[n] ;
	}
	return *this ;
}

// subtraction of two vectors
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator-( const Vector<T,N> & y )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = x[n] - y[n] ;
	}
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator-=( const Vector<T,N> & y )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] -= y[n] ;
	}
	return *this ;
}

// elementwise multiplication of two vectors
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator*( const Vector<T,N> & y )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = x[n] * y[n] ;
	}
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator*=( const Vector<T,N> & y )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] *= y[n] ;
	}
	return *this ;
}

// elementwise division of two vectors
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator/( const Vector<T,N> & y )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = x[n] / y[n] ;
	}
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator/=( const Vector<T,N> & y )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] /= y[n] ;
	}
	return *this ;
}

// scalar addition to vector
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator+( const T a )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = x[n] + a ;
	}
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator+=( const T a )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] += a ;
	}
	return *this ;
}

// scalar subtraction from vector
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator-( const T a )
const
{
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = x[n] - a ;
	}
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator-=( const T a )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] -= a ;
	}
	return *this ;
}

// scalar multiplication of vector
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator*( const T a )
const
{
	Vector<T,N> v( *this ) ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] *= a ;
	}
	/*
	Vector<T,N> v ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		v[n] = x[n] * a ;
	}
	//*/
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator*=( const T a )
{
	for( std::size_t n = 0 ; n < N ; ++n ) {
		x[n] *= a ;
	}
	return *this ;
}

// scalar division of vector
template< class T , std::size_t N >
inline Vector<T,N> Vector<T,N>::operator/( const T a )
const
{
	T b = static_cast<T>(1) / a ;
	Vector<T,N> v = (*this) * b ;
	return v ;
}
template< class T , std::size_t N >
inline Vector<T,N> & Vector<T,N>::operator/=( const T a )
{
	T b = static_cast<T>(1) / a ;
	(*this) *= b ;
	return *this ;
}

/// EXPLICIT TEMPLATE INSTANTIATIONS ///
template class Vector<float ,2> ;
template class Vector<float ,3> ;
template class Vector<float ,4> ;
template class Vector<double,2> ;
template class Vector<double,3> ;
template class Vector<double,4> ;
template class Vector<bool  ,2> ;
template class Vector<bool  ,3> ;
template class Vector<bool  ,4> ;
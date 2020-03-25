//
//  klein_summer.hpp
//  springs
//
//  Created by Jacob Herrmann on 22/03/2020.
//  Copyright © 2020 Jacob Herrmann. All rights reserved.
//

#ifndef klein_summer_hpp
#define klein_summer_hpp

template< class T >
class KleinSummer {
private:
	T s   = static_cast<T>(0.0) ;
	T cs  = static_cast<T>(0.0) ;
	T ccs = static_cast<T>(0.0) ;

public:
	KleinSummer() ;
	KleinSummer( const T ) ;
	void reset( void ) ;
	void reset( const T ) ;
	void add( const T ) ;
	T result( void ) ;
} ;

#endif /* klein_summer_hpp */


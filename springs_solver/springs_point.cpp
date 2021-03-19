#include "springs.hpp"

template< class T , std::size_t N >
void Point<T,N>::move( const Vector<T,N> & displacement )
{
	if( !this->fixed_all_dim ) {
		this->position += displacement.zero_where( this->fixed_dim ) ;
	}
}

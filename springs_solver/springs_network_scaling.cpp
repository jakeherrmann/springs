#include "springs.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  ██████  ███████ ████████         ███████  ██████  █████  ██      ███████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██ 
// ██       ██         ██            ██      ██      ██   ██ ██      ██              ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██  
// ██   ███ █████      ██            ███████ ██      ███████ ██      █████           ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████   
// ██    ██ ██         ██                 ██ ██      ██   ██ ██      ██              ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██  
//  ██████  ███████    ██    ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::get_scale_network( void )
{
	// want to reset minimum positions to the origin
	scale_position_min = nodes.begin()->point->position ;
	scale_position_max = nodes.begin()->point->position ;
	iterNode n ;
	iterNode n_end ;
	for( n = nodes.begin() , n_end = nodes.end() ; n != n_end ; ++n ) {
		for( std::size_t d = 0 ; d < N ; ++d ) {
			if( n->point->position[d] < scale_position_min[d] ) {
				scale_position_min[d] = n->point->position[d] ;
			} else if( n->point->position[d] > scale_position_max[d] ) {
				scale_position_max[d] = n->point->position[d] ;
			}
		}
	}
	scale_position_range = scale_position_max - scale_position_min ;

	// want to scale network size and stiffnesses by average spring length and stiffness
	scale_stiffness = static_cast<T>(0.0) ;
	scale_length = static_cast<T>(0.0) ;
	iterSpring s ;
	iterSpring s_end ;
	for( s = springs.begin() , s_end = springs.end() ; s != s_end ; ++s ) {
		s->spring_energy() ;
		scale_stiffness += s->spring_get_scale_stiffness() ;
		scale_length += s->length ;
	}
	scale_stiffness /= static_cast<T>(springs.size()) ;
	scale_length    /= static_cast<T>(springs.size()) ;
	scale_stiffness = ( scale_stiffness > 0.0 ) ? scale_stiffness : static_cast<T>(1.0) ;
	scale_length    = ( scale_length    > 0.0 ) ? scale_length    : static_cast<T>(1.0) ;
	scale_force = scale_length * scale_stiffness ;
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ██████  ███████ ███████  ██████  █████  ██      ███████         ███    ██ ███████ ████████ ██     ██  ██████  ██████  ██   ██ 
// ██   ██ ██      ██      ██      ██   ██ ██      ██              ████   ██ ██         ██    ██     ██ ██    ██ ██   ██ ██  ██  
// ██████  █████   ███████ ██      ███████ ██      █████           ██ ██  ██ █████      ██    ██  █  ██ ██    ██ ██████  █████   
// ██   ██ ██           ██ ██      ██   ██ ██      ██              ██  ██ ██ ██         ██    ██ ███ ██ ██    ██ ██   ██ ██  ██  
// ██   ██ ███████ ███████  ██████ ██   ██ ███████ ███████ ███████ ██   ████ ███████    ██     ███ ███   ██████  ██   ██ ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::rescale_network( const int & direction )
{
	T scale_length_inv ;
	T scale_force_inv ;
	T scale_stiffness_inv ;
	switch( direction ) {
		case +1:
			// set to normalize scale from original scale
			scale_length_inv = static_cast<T>(1) / scale_length ;
			scale_force_inv = static_cast<T>(1) / scale_force ;
			scale_stiffness_inv = static_cast<T>(1) / scale_stiffness ;
			for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
				p->position -= scale_position_min ;
				p->position *= scale_length_inv ;
				p->force_applied *= scale_force_inv ;
			}
			for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
				s->spring_rescale( scale_length_inv , scale_stiffness_inv ) ;
				//s->precompute_parameters() ;
			}
			break ;

		case -1:
			// reset to original scale from normalize scale
			for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
				p->position *= scale_length ;
				p->position += scale_position_min ;
				p->force_applied *= scale_force ;
			}
			for( iterSpring s = springs.begin() ; s != springs.end() ; ++s ) {
				s->spring_rescale( scale_length , scale_stiffness ) ;
				//s->precompute_parameters() ;
			}
			break ;

		default:
			break ;
	}
	return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ███████ ██ ███    ██ ██████          ███    ███  █████  ██   ██         ███████ ██████  ██████  ██ ███    ██  ██████          ██      ███████ ███    ██  ██████  ████████ ██   ██ 
// ██      ██ ████   ██ ██   ██         ████  ████ ██   ██  ██ ██          ██      ██   ██ ██   ██ ██ ████   ██ ██               ██      ██      ████   ██ ██          ██    ██   ██ 
// █████   ██ ██ ██  ██ ██   ██         ██ ████ ██ ███████   ███           ███████ ██████  ██████  ██ ██ ██  ██ ██   ███         ██      █████   ██ ██  ██ ██   ███    ██    ███████ 
// ██      ██ ██  ██ ██ ██   ██         ██  ██  ██ ██   ██  ██ ██               ██ ██      ██   ██ ██ ██  ██ ██ ██    ██         ██      ██      ██  ██ ██ ██    ██    ██    ██   ██ 
// ██      ██ ██   ████ ██████  ███████ ██      ██ ██   ██ ██   ██ ███████ ███████ ██      ██   ██ ██ ██   ████  ██████  ███████ ███████ ███████ ██   ████  ██████     ██    ██   ██ 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< class T , std::size_t N >
void SpringNetwork<T,N>::find_max_spring_length( void )
{
	Vector<T,N> x_min = points.begin()->position ;
	Vector<T,N> x_max = points.begin()->position ;
	for( iterPoint p = points.begin() ; p != points.end() ; ++p ) {
		for( std::size_t n = 0 ; n < N ; ++n ) {
			x_min[n] = ( p->position[n] < x_min[n] ) ? p->position[n] : x_min[n] ;
			x_max[n] = ( p->position[n] > x_max[n] ) ? p->position[n] : x_max[n] ;
		}
	}
	max_spring_length = static_cast<T>(0) ;
	T dx ;
	for( std::size_t n = 0 ; n < N ; ++n ) {
		dx = x_max[n] - x_min[n] ;
		max_spring_length = ( dx > max_spring_length ) ? dx : max_spring_length ;
	}
	max_spring_length *= 3.0 ;
	return ;
}
#include "springs.hpp"
#include <string>
#include <sys/stat.h>

void make_dir( const std::string & dir_name )
{
	struct stat info ;
	if( stat( dir_name.c_str() , &info ) != 0 ) {
		std::system( ("mkdir " + dir_name).c_str() ) ;
	} else if( info.st_mode & S_IFDIR ) {
		// directory already exists
	} else {
		// exists but not a directory
	}
	return ;
}
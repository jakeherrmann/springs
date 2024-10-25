#include "springs.hpp"
#include <string>
#include <sys/stat.h>
#include <stdexcept>

#ifdef _WIN32
	#define FILESEP '\\'
#else
	#define FILESEP '/'
#endif

void make_dir( const std::string & dir_name )
{
	struct stat info ;
	if( stat( dir_name.c_str() , &info ) != 0 ) {
		// directory does not exist --> make new
		std::system( ("mkdir " + dir_name).c_str() ) ;
	} else if( info.st_mode & S_IFDIR ) {
		// directory already exists
	} else {
		// exists but not a directory
		throw std::runtime_error( "Error: Specified output directory exists but is not a directory: " + dir_name ) ;
	}
	return ;
}
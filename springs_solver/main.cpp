//
//  main.cpp
//  springs
//
//  Created by Jacob Herrmann on 12/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <string>
#include <memory>

#ifdef _WIN32
	#include <direct.h>
	#define get_current_dir _getcwd
#else
	#include <unistd.h>
	#define get_current_dir getcwd
#endif

#ifdef _WIN32
	#define FILESEP '\\'
#else
	#define FILESEP '/'
#endif

#include "springs.hpp"

int main( int argc , const char * argv[] ) {
	// REQUIRED ARGUMENTS
	// argv[0] := PATH TO THIS SOLVER EXECUTABLE
	// argv[1] := INPUT DIRECTORY
	// argv[2] := OUTPUT DIRECTORY

	// DEFAULT PARAMETERS
	bool verbose = false ;
	bool parallel = true ;

	// OPTIONAL ARGUMENTS (DEFAULT)
	// -v --verbose := (0)=quiet, 1=verbose
	if( (argc<2) || (argc%2!=1) ) {
		std::cout << "ERROR: Incorrect number of parameter-value pairs in argument list.  Exiting." << std::endl ;
		return -1 ;
	}
	std::string arg ;
	std::string val ;
	for( int i = 3 ; i < argc ; i+=2 ) {
		arg = argv[i  ] ;
		val = argv[i+1] ;
		if( (arg=="-v") || (arg=="--verbose") ) {
			verbose = std::stoi(val) != 0 ;
		} else if( (arg=="-p") || (arg=="--parallel") ) {
			parallel = std::stoi(val) != 0 ;
		}
	}

	// if quiet mode, silence the standard output
	if( !verbose ) {
		std::cout.setstate( std::ios_base::badbit ) ;
	} else {
		std::cout.clear() ;
	}

	/*
	char cwd [FILENAME_MAX] ;
	if( get_current_dir(cwd,sizeof(cwd)) == NULL ) {
		std::cout << "COULD NOT GET CURRENT DIRECTORY" << std::endl ;
	} else {
		std::cout << "WORKING DIRECTORY:\n" << cwd << std::endl ;
	}
	//*/

	// input/output locations
	std::string dir_input  = std::string("")+".."+FILESEP+".."+FILESEP+"INPUT"+FILESEP ;
	std::string dir_output ;
	if( argc > 1 ) {
		dir_input = std::string( argv[1] ) ;
		if( dir_input.back() != FILESEP ) {
			dir_input += FILESEP ;
		}
	}
	if( argc > 2 ) {
		dir_output = std::string( argv[2] ) ;
		if( dir_output.back() != FILESEP ) {
			dir_output += FILESEP ;
		}
	} else {
		dir_output = dir_input.substr(0,dir_input.length()-1) + "_OUTPUT"+FILESEP ;
	}
	
	std::cout << "\n0.  Parameters" << std::endl ;
	NetworkParameters network_parameters( dir_input , dir_output ) ;
	network_parameters.parallelism_enabled = parallel ;
	
	std::cout << "\n1.  Create Spring Network" << std::endl ;
	std::unique_ptr<ASpringNetwork> spring_network = ASpringNetwork::create_spring_network_obj( network_parameters ) ;
	if( spring_network == NULL ) {
		return -1 ;
	}
	
	std::cout << "\n2.  Setup" << std::endl ;
	spring_network->setup( network_parameters ) ;
	
	std::cout << "\n3.  Solve" << std::endl ;
	spring_network->solve() ;
	
	std::cout << "\n4.  Complete" << std::endl ;
    return 0 ;
}

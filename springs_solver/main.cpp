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

#ifdef WINDOWS
	#include <direct.h>
	#define get_current_dir _getcwd
#else
	#include <unistd.h>
	#define get_current_dir getcwd
#endif

#include "springs.hpp"

int main( int argc , const char * argv[] ) {
	// ARGV[0] := PATH TO THIS SOLVER EXECUTABLE
	// ARGV[1] := INPUT DIRECTORY
	// ARGV[2] := OUTPUT DIRECTORY
	// ARGV[3] := NUMBER OF ITERATIONS TO SAVE
	// ARGV[4] := 0=QUIET, 1=VERBOSE

	// assume quiet mode, turn on verbose if desired
	std::cout.setstate( std::ios_base::badbit ) ;
	if( argc > 4 ) {
		if( !( argv[4][0] != '0' ) ) {
			std::cout.clear() ;
		}
	}

	// number of intermediate iterations to save
	std::size_t num_iter_save = 0 ;
	if( argc > 3 ) {
		num_iter_save = static_cast<std::size_t>( std::atoi(argv[3]) ) ;
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
	std::string dir_input  = "../../INPUT/" ;
	std::string dir_output ;
	if( argc > 1 ) {
		dir_input = std::string( argv[1] ) ;
		if( dir_input.back() != '/' ) {
			dir_input += '/' ;
		}
	}
	if( argc > 2 ) {
		dir_output = std::string( argv[2] ) ;
		if( dir_output.back() != '/' ) {
			dir_output += '/' ;
		}
	} else {
		dir_output = dir_input.substr(0,dir_input.length()-1) + "_OUTPUT/" ;
	}
	
	std::cout << "\n0.  Parameters" << std::endl ;
	NetworkParameters network_parameters( dir_input , dir_output ) ;
	
	std::cout << "\n1.  Create Spring Network" << std::endl ;
	std::unique_ptr<ASpringNetwork> spring_network = ASpringNetwork::create_spring_network_obj( network_parameters ) ;
	if( spring_network == NULL ) {
		return -1 ;
	}
	
	std::cout << "\n2.  Setup" << std::endl ;
	spring_network->setup( network_parameters ) ;
	
	std::cout << "\n3.  Solve" << std::endl ;
	spring_network->num_iter_save = num_iter_save ;
	spring_network->solve() ;
	
	std::cout << "\n4.  Complete" << std::endl ;
    return 0 ;
}

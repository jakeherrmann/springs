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
	
	char cwd [FILENAME_MAX] ;
	if( get_current_dir(cwd,sizeof(cwd)) == NULL ) {
		std::cout << "COULD NOT GET CURRENT DIRECTORY" << std::endl ;
	} else {
		std::cout << "WORKING DIRECTORY:\n" << cwd << std::endl ;
	}
	
	std::string dir_input  = "./INPUT/" ;
	std::string dir_output = "./OUTPUT/" ;
	std::system( ("mkdir " + dir_output).c_str() ) ;
	
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
	spring_network->solve() ;
	
	std::cout << "\n4.  Complete" << std::endl ;
    return 0 ;
}

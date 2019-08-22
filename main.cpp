//
//  main.cpp
//  springs
//
//  Created by Jacob Herrmann on 12/08/2019.
//  Copyright © 2019 Jacob Herrmann. All rights reserved.
//

#include <iostream>
#include <string>
#include <memory>

#include "springs.hpp"

int main( int argc , const char * argv[] ) {
	
	std::string file_parameters    = "/Users/jake/Documents/GitHub/springNet/springs/springs/INPUT/network_parameters.txt" ;
	std::string file_setup_nodes   = "/Users/jake/Documents/GitHub/springNet/springs/springs/INPUT/network_setup_nodes.dat" ;
	std::string file_setup_springs = "/Users/jake/Documents/GitHub/springNet/springs/springs/INPUT/network_setup_springs.dat" ;
	
	std::cout << "\n0.  Parameters" << std::endl ;
	NetworkParameters network_parameters( file_parameters.c_str() ) ;
	
	std::cout << "\n1.  Create Spring Network" << std::endl ;
	std::unique_ptr<ASpringNetwork> spring_network = ASpringNetwork::create_spring_network_obj( network_parameters ) ;
	if( spring_network == NULL ) {
		return -1 ;
	}
	
	std::cout << "\n2.  Setup" << std::endl ;
	spring_network->setup( network_parameters ,
						   file_setup_nodes.c_str() ,
						   file_setup_springs.c_str() ) ;
	
	std::cout << "\n3.  Execute Jobs" << std::endl ;
	spring_network->execute_jobs() ;
	
	std::cout << "Complete" << std::endl ;
	
	// read all parameters in high-memory mode
	// afterwards use if(lowmem/highmem)
	//     typecast to lower memory if option specified.
	//     then call solver of the right datatype <type>
	
    return 0 ;
}

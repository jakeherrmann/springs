#include "springs.hpp"
#include <fstream>
#include <iostream>
#include <string>

NetworkParameters::NetworkParameters( const std::string & dir_input ,
									  const std::string & dir_output )
{
	this->dir_input        = dir_input ;
	this->dir_output       = dir_output ;
	make_dir( dir_output ) ;
	file_input_parameters  = dir_input  + "network_parameters.txt" ;
	file_output_parameters = dir_output + "network_parameters.txt" ;
	load_parameters( file_input_parameters.c_str() ) ;
	save_parameters( file_output_parameters.c_str() ) ;
	return ;
}

void NetworkParameters::load_parameters( const char * file_name )
{
	std::ifstream file ;
	std::string arg ;
	std::string val ;
	if( file_name != NULL ) {
		file.open( file_name ) ;
		if( !file ) {
			std::cout << "COULD NOT OPEN FILE:\n" << file_name << std::endl ;
		} else {
			while( file >> arg >> val ) {
				if( arg=="num_points"                 ) { num_points                 = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_springs"                ) { num_springs                = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="precision"                  ) { precision                  =                                    val   ; } else
				if( arg=="num_dimensions"             ) { num_dimensions             = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="algorithm"                  ) { algorithm                  =                                    val   ; } else
				if( arg=="objective"                  ) { objective                  =                                    val   ; } else
				if( arg=="num_threads"                ) { user_num_threads           = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_save"              ) { num_iter_save              = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_print"             ) { num_iter_print             = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="num_iter_max"               ) { num_iter_max               = static_cast<std::size_t>(std::stoi(val)) ; } else
				if( arg=="include_force_fixed_nodes"  ) { include_force_fixed_nodes  =        static_cast<bool>(std::stod(val)) ; } else
				if( arg=="use_numerical_hessian"      ) { use_numerical_hessian      =        static_cast<bool>(std::stod(val)) ; } else
				if( arg=="tolerance_change_objective" ) { tolerance_change_objective =                          std::stod(val)  ; } else
				if( arg=="tolerance_sum_net_force"    ) { tolerance_sum_net_force    =                          std::stod(val)  ; }
				else { std::cout << "Parameter not recognized: " << arg << std::endl ; }
			}
			file.close() ;
		}
	}
	return ;
}

void NetworkParameters::save_parameters( const char * file_name )
{
	std::ofstream file ;
	if( file_name != NULL ) {
		file.open( file_name ) ;
		if( !file ) {
			std::cout << "COULD NOT OPEN FILE:\n" << file_name << std::endl ;
		} else {
			file << "num_points"                 << ' ' << num_points                 << '\n'
				 << "num_springs"                << ' ' << num_springs                << '\n'
				 << "precision"                  << ' ' << precision                  << '\n'
				 << "num_dimensions"             << ' ' << num_dimensions             << '\n'
				 << "algorithm"                  << ' ' << algorithm                  << '\n'
				 << "objective"                  << ' ' << objective                  << '\n'
				 << "num_threads"                << ' ' << user_num_threads           << '\n'
				 << "num_iter_save"              << ' ' << num_iter_save              << '\n'
				 << "num_iter_print"             << ' ' << num_iter_print             << '\n'
				 << "num_iter_max"               << ' ' << num_iter_max               << '\n'
				 << "include_force_fixed_nodes"  << ' ' << include_force_fixed_nodes  << '\n'
				 << "use_numerical_hessian"      << ' ' << use_numerical_hessian      << '\n'
				 << "tolerance_change_objective" << ' ' << tolerance_change_objective << '\n'
				 << "tolerance_sum_net_force"    << ' ' << tolerance_sum_net_force ;
			file.close() ;
		}
	}
	return ;
}
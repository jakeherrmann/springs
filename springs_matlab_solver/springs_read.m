function [ nodes , springs , network_param ] = springs_read( dir_output )

file_parameters = fullfile( dir_output , 'network_parameters.txt' ) ;
file_nodes      = fullfile( dir_output , 'network_nodes.dat' ) ;
file_springs    = fullfile( dir_output , 'network_springs.dat' ) ;

numeric_args = {
	'num_points'
	'num_springs'
	'num_dimensions'
	'num_threads'
	'num_stiffness_tension'
	'num_stiffness_compression'
	'num_iter_save'
	'num_iter_print'
	'num_iter_max'
	'tolerance_change_objective'
	'tolerance_sum_net_force'
	} ;
bool_args = {
	'include_force_fixed_nodes'
	'use_numerical_hessian'
	} ;
	
filename = file_parameters ;
fid = fopen( filename ,'rt') ;
while ~feof( fid )
	[ arg , val ] = strtok( fgetl( fid ) ) ;
	val = strtrim(val) ;
	if ismember( arg , numeric_args )
		val = str2num(val) ;
	elseif ismember( arg , bool_args )
		val = logical(str2num(val)) ;
	end
	network_param.(arg) = val ; 
end
fclose( fid ) ;

nodes.position = zeros( [ network_param.num_points , network_param.num_dimensions ] ) ;
nodes.force    = zeros( [ network_param.num_points , network_param.num_dimensions ] ) ;
nodes.fixed    = false( [ network_param.num_points , network_param.num_dimensions ] ) ;

springs.nodes                               = zeros( [ network_param.num_springs , 2 ] ) ;
springs.rest_length                         = zeros( [ network_param.num_springs , 1 ] ) ;
springs.force_length_type_tension           = zeros( [ network_param.num_springs , 1 ] ) ;
springs.force_length_type_compression       = zeros( [ network_param.num_springs , 1 ] ) ;
springs.force_length_parameters_tension     =  cell( [ network_param.num_springs , 1 ] ) ;
springs.force_length_parameters_compression =  cell( [ network_param.num_springs , 1 ] ) ;

filename = file_nodes ;
fid = fopen( filename ,'rb') ;
for pp = 1 : network_param.num_points
	nodes.position(pp,:) = fread( fid , network_param.num_dimensions , network_param.precision ) ;
	nodes.force(pp,:)    = fread( fid , network_param.num_dimensions , network_param.precision ) ;
	nodes.fixed(pp,:)    = fread( fid , network_param.num_dimensions , 'uint8'                 ) ;
end
fclose( fid ) ;

filename = file_springs ;
fid = fopen( filename ,'rb') ;
for ss = 1 : network_param.num_springs
	springs.nodes(ss,:)                             = fread( fid , 2     , 'uint32'                ) + 1 ;
	springs.force_length_type_tension(ss)           = fread( fid , 1     , 'uint8'                 ) ;
	springs.force_length_type_compression(ss)       = fread( fid , 1     , 'uint8'                 ) ;
	springs.rest_length(ss)                         = fread( fid , 1     , network_param.precision ) ;
	NFLPT                                           = fread( fid , 1     , 'uint32'                ) ;
	springs.force_length_parameters_tension{ss}     = fread( fid , NFLPT , network_param.precision ) ;
	NFLPC                                           = fread( fid , 1     , 'uint32'                ) ;
	springs.force_length_parameters_compression{ss} = fread( fid , NFLPC , network_param.precision ) ;
end
fclose( fid ) ;

end

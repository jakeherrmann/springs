function [ nodes , springs , network_param ] = springs_read( dir_output )

file_parameters = fullfile( dir_output , 'network_parameters.txt' ) ;
file_nodes      = fullfile( dir_output , 'network_nodes.dat' ) ;
file_springs    = fullfile( dir_output , 'network_springs.dat' ) ;

numeric_args = {
	'num_points'
	'num_springs'
	'num_dimensions'
	'num_stiffness_tension'
	'num_stiffness_compression'
	'num_iter_save'
	'num_iter_print'
	'num_iter_max'
	'tolerance_change_energy'
	'tolerance_sum_net_force'
	} ;
	
filename = file_parameters ;
fid = fopen( filename ,'rt') ;
while ~feof( fid )
	[ arg , val ] = strtok( fgetl( fid ) ) ;
	val = strtrim(val) ;
	if ismember( arg , numeric_args )
		val = str2num(val) ;
	end
	network_param.(arg) = val ; 
end
fclose( fid ) ;

nodes.position = zeros( [ network_param.num_points , network_param.num_dimensions ] ) ;
nodes.force    = zeros( [ network_param.num_points , network_param.num_dimensions ] ) ;
nodes.fixed    = false( [ network_param.num_points , 1                            ] ) ;

springs.nodes                 = zeros( [ network_param.num_springs , 2                                       ] ) ;
springs.rest_length           = zeros( [ network_param.num_springs , 1                                       ] ) ;
springs.stiffness_tension     = zeros( [ network_param.num_springs , network_param.num_stiffness_tension     ] ) ;
springs.stiffness_compression = zeros( [ network_param.num_springs , network_param.num_stiffness_compression ] ) ;
springs.compression           = zeros( [ network_param.num_springs , 1                                       ] ) ;

filename = file_nodes ;
fid = fopen( filename ,'rb') ;
for pp = 1 : network_param.num_points
	nodes.position(pp,:) = fread( fid , network_param.num_dimensions , network_param.precision ) ;
	nodes.force(pp,:)    = fread( fid , network_param.num_dimensions , network_param.precision ) ;
	nodes.fixed(pp)      = fread( fid , 1                            , 'uint8'                 ) ;
end
fclose( fid ) ;

filename = file_springs ;
fid = fopen( filename ,'rb') ;
for ss = 1 : network_param.num_springs
	springs.nodes(ss,:)                 = fread( fid , 2                                       , 'uint32'                ) + 1 ;
	springs.rest_length(ss)             = fread( fid , 1                                       , network_param.precision ) ;
	springs.stiffness_tension(ss,:)     = fread( fid , network_param.num_stiffness_tension     , network_param.precision ) ;
	springs.stiffness_compression(ss,:) = fread( fid , network_param.num_stiffness_compression , network_param.precision ) ;
	springs.compression(ss)             = fread( fid , 1                                       , 'uint8'                 ) ;
end
fclose( fid ) ;

end

function [ network_parameters , nodes , springs ] = read_spring_network( dir_output )

file_parameters = fullfile( dir_output , 'network_parameters.txt' ) ;
file_nodes      = fullfile( dir_output , 'network_nodes.dat' ) ;
file_springs    = fullfile( dir_output , 'network_springs.dat' ) ;

filename = file_parameters ;
fid = fopen( filename ,'rt') ;
num_points                = str2num( fgetl( fid ) ) ;
num_springs               = str2num( fgetl( fid ) ) ;
precision                 =          fgetl( fid )   ;
num_dimensions            = str2num( fgetl( fid ) ) ;
num_stiffness_tension     = str2num( fgetl( fid ) ) ;
num_stiffness_compression = str2num( fgetl( fid ) ) ;
fclose( fid ) ;

nodes.position   = zeros( [ num_points , num_dimensions ] ) ;
nodes.force      = zeros( [ num_points , num_dimensions ] ) ;
nodes.fixed      = false( [ num_points , 1              ] ) ;
nodes.referenced = false( [ num_points , 1              ] ) ;
fid = fopen( file_nodes ,'rb') ;
for pp = 1 : num_points
	nodes.position(pp,:) = fread( fid , num_dimensions , precision ) ;
	nodes.force(pp,:)    = fread( fid , num_dimensions , precision ) ;
	nodes.fixed(pp)      = fread( fid , 1              , 'uint8'   ) ;
	nodes.referenced(pp) = fread( fid , 1              , 'uint8'   ) ;
end
fclose( fid ) ;

springs.nodes                 = zeros( [ num_springs , 2                         ] ) ;
springs.rest_length           = zeros( [ num_springs , 1                         ] ) ;
springs.stiffness_tension     = zeros( [ num_springs , num_stiffness_tension     ] ) ;
springs.stiffness_compression = zeros( [ num_springs , num_stiffness_compression ] ) ;
springs.compressible          = zeros( [ num_springs , 1                         ] ) ;
springs.force                 = zeros( [ num_springs , 1                         ] ) ;
springs.strain                = zeros( [ num_springs , 1                         ] ) ;
springs.broken                = false( [ num_springs , 1                         ] ) ;
springs.repairable            = false( [ num_springs , 1                         ] ) ;
fid = fopen( file_springs ,'rb') ;
for ss = 1 : num_springs
	springs.nodes(ss,:)                 = fread( fid , 2                         , 'uint32'  ) + 1 ;
	springs.rest_length(ss)             = fread( fid , 1                         , precision ) ;
	springs.stiffness_tension(ss,:)     = fread( fid , num_stiffness_tension     , precision ) ;
	springs.stiffness_compression(ss,:) = fread( fid , num_stiffness_compression , precision ) ;
	springs.compressible(ss)            = fread( fid , 1                         , 'uint8'   ) ;
	springs.force(ss)                   = fread( fid , 1                         , precision ) ;
	springs.strain(ss)                  = fread( fid , 1                         , precision ) ;
	springs.broken(ss)                  = fread( fid , 1                         , 'uint8'   ) ;
	springs.repairable(ss)              = fread( fid , 1                         , 'uint8'   ) ;
end
fclose( fid ) ;

network_parameters = struct( ...
	'num_points'                , num_points                , ...
	'num_springs'               , num_springs               , ...
	'precision'                 , precision                 , ...
	'num_dimensions'            , num_dimensions            , ...
	'num_stiffness_tension'     , num_stiffness_tension     , ...
	'num_stiffness_compression' , num_stiffness_compression ) ;

end

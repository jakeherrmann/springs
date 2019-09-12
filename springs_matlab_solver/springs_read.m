function [ nodes , springs ] = springs_read( dir_output )

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

nodes.position = zeros( [ num_points , num_dimensions ] ) ;
nodes.force    = zeros( [ num_points , num_dimensions ] ) ;
nodes.fixed    = false( [ num_points , 1              ] ) ;

springs.nodes                 = zeros( [ num_springs , 2                         ] ) ;
springs.rest_length           = zeros( [ num_springs , 1                         ] ) ;
springs.stiffness_tension     = zeros( [ num_springs , num_stiffness_tension     ] ) ;
springs.stiffness_compression = zeros( [ num_springs , num_stiffness_compression ] ) ;
springs.compression           = zeros( [ num_springs , 1                         ] ) ;

filename = file_nodes ;
fid = fopen( filename ,'rb') ;
for pp = 1 : num_points
	nodes.position(pp,:) = fread( fid , num_dimensions , precision ) ;
	nodes.force(pp,:)    = fread( fid , num_dimensions , precision ) ;
	nodes.fixed(pp)      = fread( fid , 1              , 'uint8'   ) ;
end
fclose( fid ) ;

filename = file_springs ;
fid = fopen( filename ,'rb') ;
for ss = 1 : num_springs
	springs.nodes(ss,:)                 = fread( fid , 2                         , 'uint32'  ) + 1 ;
	springs.rest_length(ss)             = fread( fid , 1                         , precision ) ;
	springs.stiffness_tension(ss,:)     = fread( fid , num_stiffness_tension     , precision ) ;
	springs.stiffness_compression(ss,:) = fread( fid , num_stiffness_compression , precision ) ;
	springs.compression(ss)             = fread( fid , 1                         , 'uint8'   ) ;
end
fclose( fid ) ;

end

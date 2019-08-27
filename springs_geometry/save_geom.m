function save_geom( dir_input , precision , nodes , springs )

if ~exist( dir_input ,'dir')
	mkdir( dir_input ) ;
end

num_points     = size( nodes.position ,1) ;
num_dimensions = size( nodes.position ,2) ;
num_springs    = size( springs.nodes   ,1) ;
num_stiffness_tension     = size( springs.stiffness_tension     ,2) ;
num_stiffness_compression = size( springs.stiffness_compression ,2) ;

filename = fullfile( dir_input , 'network_parameters.txt' ) ;
fid = fopen( filename ,'wt') ;
fprintf( fid , '%d\n' , num_points ) ;
fprintf( fid , '%d\n' , num_springs ) ;
fprintf( fid , '%s\n' , precision ) ;
fprintf( fid , '%d\n' , num_dimensions ) ;
fprintf( fid , '%d\n' , num_stiffness_tension ) ;
fprintf( fid , '%d\n' , num_stiffness_compression ) ;
fclose( fid ) ;

filename = fullfile( dir_input , 'network_nodes.dat' ) ;
fid = fopen( filename ,'wb') ;
for nn = 1 : num_points
	fwrite( fid , nodes.position(nn,:) , precision ) ;
	fwrite( fid , nodes.force(nn,:)    , precision ) ;
	fwrite( fid , nodes.fixed(nn)      , 'uint8'   ) ;
end
fclose( fid ) ;

filename = fullfile( dir_input , 'network_springs.dat' ) ;
fid = fopen( filename ,'wb') ;
for ss = 1 : num_springs
	fwrite( fid , springs.nodes(ss,:)-1               , 'uint32'  ) ;
	fwrite( fid , springs.rest_length(ss)             , precision ) ;
	fwrite( fid , springs.stiffness_tension(ss,:)     , precision ) ;
	fwrite( fid , springs.stiffness_compression(ss,:) , precision ) ;
	fwrite( fid , springs.compression(ss)             , 'uint8'   ) ;
end
fclose( fid ) ;

end
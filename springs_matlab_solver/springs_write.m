function springs_write( dir_input , nodes , springs , varargin )

if ~exist( dir_input ,'dir')
	mkdir( dir_input ) ;
end

if numel(varargin) > 0
	for ff = fieldnames(varargin{1})'
		network_param.(ff{1}) = varargin{1}.(ff{1}) ;
	end
end

network_param.num_points                = size( nodes.position ,1) ;
network_param.num_springs               = size( springs.nodes  ,1) ;
network_param.num_dimensions            = size( nodes.position ,2) ;
network_param.num_stiffness_tension     = size( springs.stiffness_tension     ,2) ;
network_param.num_stiffness_compression = size( springs.stiffness_compression ,2) ;

numeric_args = {
	'num_points'
	'num_springs'
	'num_dimensions'
	'num_stiffness_tension'
	'num_stiffness_compression'
	'num_iter_save'
	'num_iter_print'
	'num_iter_max'
	'use_sum_net_force'
	'use_numerical_hessian'
	'tolerance_change_energy'
	'tolerance_sum_net_force'
	} ;

filename = fullfile( dir_input , 'network_parameters.txt' ) ;
fid = fopen( filename ,'wt') ;
for ff = fieldnames(network_param)'
	if ismember( ff{1} , numeric_args )
		fprintf( fid , '%s %s\n' , ff{1} , num2str( network_param.(ff{1}) ) ) ;
	else
		fprintf( fid , '%s %s\n' , ff{1} , network_param.(ff{1}) ) ;
	end
end
fclose( fid ) ;

filename = fullfile( dir_input , 'network_nodes.dat' ) ;
fid = fopen( filename ,'wb') ;
for nn = 1 : network_param.num_points
	fwrite( fid , nodes.position(nn,:) , network_param.precision ) ;
	fwrite( fid , nodes.force(nn,:)    , network_param.precision ) ;
	fwrite( fid , nodes.fixed(nn,:)    , 'uint8'                 ) ;
end
fclose( fid ) ;

filename = fullfile( dir_input , 'network_springs.dat' ) ;
fid = fopen( filename ,'wb') ;
for ss = 1 : network_param.num_springs
	fwrite( fid , springs.nodes(ss,:)-1               , 'uint32'                ) ;
	fwrite( fid , springs.rest_length(ss)             , network_param.precision ) ;
	fwrite( fid , springs.stiffness_tension(ss,:)     , network_param.precision ) ;
	fwrite( fid , springs.stiffness_compression(ss,:) , network_param.precision ) ;
	fwrite( fid , springs.compression(ss)             , 'uint8'                 ) ;
end
fclose( fid ) ;

end
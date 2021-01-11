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

numeric_args = {
	'num_points'
	'num_springs'
	'num_dimensions'
	'num_threads'
	'num_iter_save'
	'num_iter_print'
	'num_iter_max'
	'include_force_fixed_nodes'
	'use_numerical_hessian'
	'tolerance_change_objective'
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
	fwrite( fid , springs.nodes(ss,:)-1                                  , 'uint32'                ) ;
	fwrite( fid , springs.force_length_type_tension(ss)                  , 'uint8'                 ) ;
	fwrite( fid , springs.force_length_type_compression(ss)              , 'uint8'                 ) ;
	fwrite( fid , springs.rest_length(ss)                                , network_param.precision ) ;
	fwrite( fid , numel(springs.force_length_parameters_tension{ss}    ) , 'uint32'                ) ;
	fwrite( fid , springs.force_length_parameters_tension{ss}            , network_param.precision ) ;
	fwrite( fid , numel(springs.force_length_parameters_compression{ss}) , 'uint32'                ) ;
	fwrite( fid , springs.force_length_parameters_compression{ss}        , network_param.precision ) ;
end
fclose( fid ) ;

end
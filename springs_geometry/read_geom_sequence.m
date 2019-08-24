
clear
clc

dir_input = '/Users/jake/Documents/GitHub/springNet/springs/INPUT/' ;
dir_output = strrep( dir_input , 'INPUT' , 'OUTPUT' ) ;

file_parameters = fullfile( dir_input  , 'network_parameters.txt' ) ;
% file_nodes      = fullfile( dir_output , 'network_output_nodes.dat' ) ;
% file_springs    = fullfile( dir_output , 'network_output_springs.dat' ) ;

filename = file_parameters ;
fid = fopen( filename ,'rt') ;
num_points                = str2num( fgetl( fid ) ) ;
num_springs               = str2num( fgetl( fid ) ) ;
precision                 =          fgetl( fid )   ;
num_dimensions            = str2num( fgetl( fid ) ) ;
num_stiffness_tension     = str2num( fgetl( fid ) ) ;
num_stiffness_compression = str2num( fgetl( fid ) ) ;
fclose( fid ) ;

file_list_nodes = dir( fullfile( dir_output , 'iter_*_nodes.dat' ) ) ;
temp = str2double( strrep(strrep( { file_list_nodes.name }' ,'iter_',''),'_nodes.dat','') ) ;
[ ~ , ind_sort ] = sort( temp ) ;
file_list_nodes = file_list_nodes( ind_sort ) ;
num_file = numel( file_list_nodes ) ;
num_file = min( num_file , 100 ) ;
file_list_nodes = file_list_nodes( 1 : num_file ) ;

nodes.position = zeros( [ num_points , num_dimensions , num_file ] ) ;
nodes.force = zeros( [ num_points , num_dimensions , num_file ] ) ;
nodes.fixed = false( [ num_points , 1              , num_file ] ) ;

springs.nodes                 = zeros( [ num_springs , 2                         , num_file ] ) ;
springs.restlength            = zeros( [ num_springs , 1                         , num_file ] ) ;
springs.stiffness_tension     = zeros( [ num_springs , num_stiffness_tension     , num_file ] ) ;
springs.stiffness_compression = zeros( [ num_springs , num_stiffness_compression , num_file ] ) ;
springs.compression           = zeros( [ num_springs , 1                         , num_file ] ) ;

for ff = 1 : num_file
	filename = fullfile( dir_output , file_list_nodes(ff).name ) ;
	fid = fopen( filename ,'rb') ;
	for pp = 1 : num_points
		nodes.position(pp,:,ff) = fread( fid , num_dimensions , precision ) ;
		nodes.force(pp,:,ff) = fread( fid , num_dimensions , precision ) ;
		nodes.fixed(pp,:,ff) = fread( fid , 1              , 'uint8'   ) ;
	end
	fclose( fid ) ;
	
	filename = strrep( filename , '_nodes' , '_springs' ) ;
	fid = fopen( filename ,'rb') ;
	for ss = 1 : num_springs
		springs.nodes(ss,:,ff)                 = fread( fid , 2                         , 'uint32'  ) + 1 ;
		springs.restlength(ss,:,ff)            = fread( fid , 1                         , precision ) ;
		springs.stiffness_tension(ss,:,ff)     = fread( fid , num_stiffness_tension     , precision ) ;
		springs.stiffness_compression(ss,:,ff) = fread( fid , num_stiffness_compression , precision ) ;
		springs.compression(ss,:,ff)           = fread( fid , 1                         , 'uint8'   ) ;
	end
	fclose( fid ) ;
	
end

%%

switch num_dimensions
	case 2
		close all
		figure
		axes( ...'XLim' , [ -5 15 ]+4 , ...'YLim' , [ -5 15 ] , ...
			'DataAspectRatio' , [1,1,1] , ...
			'NextPlot' , 'add' )
		ff = 1 ;
		for ss = 1 : num_springs
			hplot_springs(ss) = plot( ...
				nodes.position( springs.nodes(ss,:,ff) ,1,ff) , ...
				nodes.position( springs.nodes(ss,:,ff) ,2,ff) , ...
				'k-' , ...
				'LineWidth' , 2 ) ;
		end
		for pp = 1 : num_points
			hplot_forces(pp) = plot( ...
				nodes.position(pp,1,ff) + [0,nodes.force(pp,1,ff)] , ...
				nodes.position(pp,2,ff) + [0,nodes.force(pp,2,ff)] , ...
				'r:' , ...
				'LineWidth' , 2 ) ;
		end
		for ff = 1 : num_file
			for ss = 1 : num_springs
				set( hplot_springs(ss) , ...
					'XData' , nodes.position( springs.nodes(ss,:,ff) ,1,ff) , ...
					'YData' , nodes.position( springs.nodes(ss,:,ff) ,2,ff) )
			end
			for pp = 1 : num_points
				set( hplot_forces(pp) , ...
					'XData' , nodes.position(pp,1,ff) + [0,nodes.force(pp,1,ff)] , ...
					'YData' , nodes.position(pp,2,ff) + [0,nodes.force(pp,2,ff)] )
			end
			title( sprintf( '%5d' , ff ) )
			pause(0.1)
			drawnow()
		end
		
	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
	case 3
		saving = false ;
		force_scaling = 1/5 ;
		close all
		figure( ...
			'Position' , [ 0 , 0 , 800 , 700 ] , ...
			'Color' , [1,1,1] ) ;
		axes( ...
			'XLim' , [ -4 , 8 ] , ...
			'YLim' , [ -4 , 8 ] , ...
			'ZLim' , [ -4 , 8 ] , ...
			'View' , [ 0 , 0 ] , ...
			'DataAspectRatio' , [1,1,1] , ...
			'Clipping' , 'off' , ...
			'Visible' , 'off' , ...
			'NextPlot' , 'add' )
		ff = 1 ;
		for ss = 1 : num_springs
			color_dim = 2 ;
			color_val = ( mean(nodes.position( springs.nodes(ss,:,ff) ,color_dim,ff)) - min(nodes.position(:,color_dim,ff)) ) / range(nodes.position(:,color_dim,ff)) ;
			color_val = color_val * 0.9 ;
			hplot_springs(ss) = plot3( ...
				nodes.position( springs.nodes(ss,:,ff) ,1,ff) , ...
				nodes.position( springs.nodes(ss,:,ff) ,2,ff) , ...
				nodes.position( springs.nodes(ss,:,ff) ,3,ff) , ...
				'k-' , ...
				'Color' , [0,0,0] + color_val , ...
				'LineWidth' , 1 ) ;
		end
		for pp = 1 : num_points
			hplot_forces(pp) = plot3( ...
				nodes.position(pp,1,ff) + [0,nodes.force(pp,1,ff)*force_scaling] , ...
				nodes.position(pp,2,ff) + [0,nodes.force(pp,2,ff)*force_scaling] , ...
				nodes.position(pp,3,ff) + [0,nodes.force(pp,3,ff)*force_scaling] , ...
				'r-' , ...
				'LineWidth' , 1 ) ;
		end
		clear frames
		for ff = 1 : num_file
			for ss = 1 : num_springs
				set( hplot_springs(ss) , ...
					'XData' , nodes.position( springs.nodes(ss,:,ff) ,1,ff) , ...
					'YData' , nodes.position( springs.nodes(ss,:,ff) ,2,ff) , ...
					'ZData' , nodes.position( springs.nodes(ss,:,ff) ,3,ff) )
			end
			for pp = 1 : num_points
				set( hplot_forces(pp) , ...
					'XData' , nodes.position(pp,1,ff) + [0,nodes.force(pp,1,ff)*force_scaling] , ...
					'YData' , nodes.position(pp,2,ff) + [0,nodes.force(pp,2,ff)*force_scaling] , ...
					'ZData' , nodes.position(pp,3,ff) + [0,nodes.force(pp,3,ff)*force_scaling] )
			end
			title( sprintf( '%5d' , ff ) )
			% set( gca , 'View' , [ 30+30*cos(2*pi*ff/num_file) , 20 ] )
			drawnow()
			if saving
				frames(ff) = getframe( gcf() ) ;
			end
		end
		
		if saving
			vid = VideoWriter( 'test.mp4' , 'MPEG-4' ) ;
			vid.FrameRate = numel(frames)/10.0 ;
			vid.open() ;
			vid.writeVideo( frames ) ;
			vid.close() ;
		end

end

%%
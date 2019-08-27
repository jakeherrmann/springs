
clear
clc

dir_output_list = dir( fullfile('..','..','TEST_RESULT','STRETCH_*') ) ;
for ii = 1 : numel( dir_output_list )
	
	%%
	dir_output = fullfile('..','..','TEST_RESULT',dir_output_list( ii ).name) ;
	
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
	
	%%
	
	switch num_dimensions
		case 2
			figure( ...
				'Position' , [0,0,1400,700] , ...
				'Color' , [1,1,1] )
			axes( ...
				'XLim' , [ -1 8 ] , ...
				'YLim' , [ -1 2 ] , ...
				'DataAspectRatio' , [1,1,1] , ...
				'TickDir' , 'out' , ...
				'NextPlot' , 'add' )
			for ss = 1 : num_springs
				plot( ...
					nodes.position( springs.nodes(ss,:) ,1) , ...
					nodes.position( springs.nodes(ss,:) ,2) , ...
					'k-' , ...
					'LineWidth' , 2 )
			end
			for pp = 1 : num_points
				plot( ...
					nodes.position(pp,1) + [0,nodes.force(pp,1)] , ...
					nodes.position(pp,2) + [0,nodes.force(pp,2)] , ...
					'r:' , ...
					'LineWidth' , 2 )
			end
			drawnow()
			frames(ii) = getframe( gcf() ) ;
			delete( gcf() )
			
		case 3
			figure( ...
				'Position' , [0,0,1400,700] , ...
				'Color' , [1,1,1] )
			axes( ...
				'XLim' , [ -1 34 ] , ...
				'YLim' , [ -1 34 ] , ...
				'ZLim' , [ -1 34 ] , ...
				'View' , [0,90],...%[ -25 , 25 ] , ...%
				'DataAspectRatio' , [1,1,1] , ...
				'NextPlot' , 'add' )
			plot3( ...
				[ nodes.position( springs.nodes(:,1) ,1) , nodes.position( springs.nodes(:,2) ,1) ]' , ...
				[ nodes.position( springs.nodes(:,1) ,2) , nodes.position( springs.nodes(:,2) ,2) ]' , ...
				[ nodes.position( springs.nodes(:,1) ,3) , nodes.position( springs.nodes(:,2) ,3) ]' , ...
				'k-' , ...
				'LineWidth' , 1 )
% 			for pp = 1 : num_points
% 				plot3( ...
% 					nodes.position(pp,1) + [0,nodes.force(pp,1)] , ...
% 					nodes.position(pp,2) + [0,nodes.force(pp,2)] , ...
% 					nodes.position(pp,3) + [0,nodes.force(pp,3)] , ...
% 					'r-' , ...
% 					'LineWidth' , 1 )
% 			end
			drawnow()
			frames(ii) = getframe( gcf() ) ;
% 			delete( gcf() )
			
	end
	
	%%
	
end

%%

return

%%

frames( arrayfun(@(f)isempty(f.cdata),frames) ) = [] ;
saving = true ;
if saving
	vid = VideoWriter( 'test.mp4' , 'MPEG-4' ) ;
	vid.FrameRate = numel(frames)/10.0 ;
	vid.open() ;
	vid.writeVideo( frames ) ;
	vid.close() ;
end
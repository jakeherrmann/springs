
dir_main = fullfile('..','..','.') ;
dir_save_list = dir(fullfile( dir_main , 'STRETCH_*' )) ;

clear frames
for dd = 1 : 130%numel(dir_save_list)
	dir_save = fullfile( dir_main , dir_save_list(dd).name ) ;
	[ nodes , springs , netparam ] = springs_read_py( dir_save ,false) ;
	switch netparam.num_dimensions
		case 2
			h = display_2D( nodes , springs , 'stiffness_tension' , [0,2] ) ;
			%ax_lims = { [-100,175] , [-100,175] , [0,1] } ;
			ax_lims = { [-50,125] , [-50,125] , [0,1] } ;
			ax_view = [0,90] ;
		case 3
			h = display_3D( nodes , springs ) ;
			ax_lims = { [-30,50] , [-30,50] , [-30,50] } ;
			ax_view = [45,45] ;
	end
	set( h.fig , ...
		'Position' , [ 0 , 0 , 750 , 750 ] )
	set( h.ax , ...
		'Position' , [ 0 , 0 , 1 , 1 ] , ...
		'Visible' , 'off' , ...
		'XLim' , ax_lims{1} , ...
		'YLim' , ax_lims{2} , ...
		'ZLim' , ax_lims{3} , ...
		'View' , ax_view )
	drawnow()
	frames(dd) = getframe( gcf() ) ;
	close( gcf() )
end
frames = frames(arrayfun(@(c)~isempty(c.cdata),frames)) ;

%%

vid = VideoWriter( fullfile('..','..','test.mp4') , 'MPEG-4' ) ;
vid.FrameRate = 10 ;
vid.open() ;
framesSubset = frames(1:2:end) ;
% framesSubset = frames(:) ;
vid.writeVideo( framesSubset ) ;
vid.close() ;
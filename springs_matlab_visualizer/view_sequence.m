
dir_main = fullfile('..','..','.') ;
dir_save_list = dir(fullfile( dir_main , 'STRETCH_*' )) ;

dir_save = fullfile( dir_main , dir_save_list(1).name ) ;
[ nodes , springs , netparam ] = springs_read_py( dir_save ,false) ;
h = display_2D( nodes , springs , 'stiffness_tension' , [0,2] ) ;
frames = getframe( gcf() ) ;
close( gcf() ) ;
for dd = 1:numel(dir_save_list)
	dir_save = fullfile( dir_main , dir_save_list(dd).name ) ;
	[ nodes , springs , netparam ] = springs_read_py( dir_save ,false) ;
	h = display_2D( nodes , springs , 'stiffness_tension' , [0,2] ) ;
	set( h.ax , ...
		'XLim' , [ -50 , 100 ] , ...
		'YLim' , [ -50 , 100 ] )
	drawnow()
	frames(dd+1) = getframe( gcf() ) ;
	close( gcf() )
end
frames = frames(arrayfun(@(c)~isempty(c.cdata),frames)) ;

%%

vid = VideoWriter( fullfile('..','..','test.mp4') , 'MPEG-4' ) ;
vid.FrameRate = 30 ;
vid.open() ;
framesSubset = frames(1:2:end) ;
% framesSubset = frames(:) ;
vid.writeVideo( framesSubset ) ;
vid.close() ;
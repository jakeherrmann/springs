
dir_main = fullfile('..','..','.') ;
dir_save_list = dir(fullfile( dir_main , 'STRETCH_*' )) ;

dir_save = fullfile( dir_main , dir_save_list(1).name ) ;
[ netparam , nodes , springs ] = read_spring_network( dir_save ) ;
h = display_2D( nodes , springs ) ;
frames = getframe( gcf() ) ;
for dd = 1:numel(dir_save_list)
	dir_save = fullfile( dir_main , dir_save_list(dd).name ) ;
	[ netparam , nodes , springs ] = read_spring_network( dir_save ) ;
	display_2D( nodes , springs , h ) ;
	drawnow()
	frames(dd+1) = getframe( gcf() ) ;
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
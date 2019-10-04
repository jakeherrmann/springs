
dir_main = fullfile('..','..') ;
dir_save_list = dir(fullfile( dir_main , 'STRETCH_*' )) ;
h = display_2D( nodes , springs ) ;
for dd = 1 : numel(dir_save_list)
	dir_save = fullfile( dir_main , dir_save_list(dd).name ) ;
	[ netparam , nodes , springs ] = read_spring_network( dir_save ) ;
	display_2D( nodes , springs , h )
	drawnow()
end
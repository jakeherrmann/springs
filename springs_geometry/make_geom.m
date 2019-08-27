
dir_input = fullfile('..','..','TEST') ;
precision = 'float' ;
geom_type = 'grid_3D' ;

switch geom_type
	case 'hexagon_2D'
		geom_size = [ 10 , 7 ] ;
		[ nodes , springs ] = make_geom_hexagon_2D( geom_size ) ;
	
	case 'random_2D'
		[ nodes , springs ] = make_geom_random_2D() ;
		
	case 'truncOct_3D'
		geom_size = [ 3 , 3 , 3 ] ;
		[ nodes , springs ] = make_geom_truncOct_3D( geom_size ) ;
		
	case 'grid_3D'
		geom_size = [10,10,10];%[ 12 , 12 , 12 ] ;
		[ nodes , springs ] = make_geom_grid_3D( geom_size ) ;
end

save_geom( dir_input , precision , nodes , springs ) ;
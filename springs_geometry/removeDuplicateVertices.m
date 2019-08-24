function B = removeDuplicateVertices( B )

[ ~ , indV , indF ] = unique( round( B.vertices ,9) ,'rows','sorted') ;
B.vertices = B.vertices(indV,:) ;
tmp = ~isnan( B.faces ) ;
B.faces(tmp) = indF( B.faces(tmp) ) ;

end
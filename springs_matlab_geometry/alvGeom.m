function [ fv ] = alvGeom()

% radius to centroid of sqr is sqrt(2/3)
% radius to centroid of hex is sqrt(5/6)

% dihedral angle between hex and hex is acos(-1/3)
% dihedral angle between hex and sqr is acos(-1/sqrt(3))

%%

% permutohedron!
n = 4 ;
pts = perms( 1:n ) ;

% unit side length
a = 1.0 ;

u = [sqrt( 2)/ 2, -sqrt( 2)/ 2,            0,          0] * a / sqrt(2) ;
v = [sqrt( 6)/ 6,  sqrt( 6)/ 6,   -sqrt(2/3),          0] * a / sqrt(2) ;
w = [sqrt(12)/12,  sqrt(12)/12,  sqrt(12)/12, -sqrt(3)/2] * a / sqrt(2) ;

x = sum( bsxfun(@times, pts , u ) ,2);
y = sum( bsxfun(@times, pts , v ) ,2);
z = sum( bsxfun(@times, pts , w ) ,2);

% vertices and faces
vertices = [ x , y , z ] ;
faces_4gon = [
	 1 ,  2 ,  8 ,  7
	 4 ,  6 , 12 , 10
	 3 ,  9 , 11 ,  5
	17 , 23 , 24 , 18
	13 , 19 , 21 , 15
	14 , 16 , 22 , 20
	] ;
faces_6gon = [
	 1 ,  3 ,  5 ,  6 ,  4 ,  2
	 1 ,  7 , 13 , 15 ,  9 ,  3
	 2 ,  4 , 10 , 16 , 14 ,  8%
	 5 , 11 , 17 , 18 , 12 ,  6%
	19 , 20 , 22 , 24 , 23 , 21%
	10 , 12 , 18 , 24 , 22 , 16%
	 9 , 15 , 21 , 23 , 17 , 11%
	 7 ,  8 , 14 , 20 , 19 , 13
	] ;

fv.sqr = struct( ...
	'faces' , faces_4gon , ...
	'vertices' , vertices ) ;
fv.hex = struct( ...
	'faces' , faces_6gon , ...
	'vertices' , vertices ) ;

%% VIEW

% patch( ...
% 	struct( ...
% 	'faces' , faces_4gon , ...
% 	'vertices' , vertices ) , ...
% 	'FaceAlpha' , 0.7 , ...
% 	'FaceColor' , [1,0.7,0.7] , ...
% 	'EdgeColor' , [0,0,0] , ...
% 	'EdgeAlpha' , 1 )
% patch( ...
% 	struct( ...
% 	'faces' , faces_6gon , ...
% 	'vertices' , vertices ) , ...
% 	'FaceAlpha' , 0.7 , ...
% 	'FaceColor' , [0.7,0.7,1] , ...
% 	'EdgeColor' , [0,0,0] , ...
% 	'EdgeAlpha' , 1 )
% text( ...
% 	vertices(:,1) , ...
% 	vertices(:,2) , ...
% 	vertices(:,3) , ...
% 	cellstr(num2str((1:24)')) , ...
% 	'FontSize' , 20 , ...
% 	'FontWeight' , 'bold' )

%% add center points to create triangles

fv.tri = struct( ...
	'faces' , [] , ...
	'vertices' , [] ) ;

for ss = {'sqr','hex'}
	for ff = 1 : size( fv.(ss{1}).faces ,1)
		face = fv.(ss{1}).faces(ff,:) ;
		centroid = mean( vertices( face ,:) ) ;
		vertices = [
			vertices
			centroid
			] ;
		cc = size(vertices,1) ;
		fv.tri.faces = [
				fv.tri.faces
				[ face.' , circshift(face.',-1) , zeros([numel(face),1])+cc ]
				] ;
	end
end
fv.tri.vertices = vertices ;

%% VIEW

% figure( ...
% 	'Color' , [1,1,1] )
% axes( ...
% 	'DataAspectRatio' , [1,1,1] , ...
% 	'Visible' , 'off' , ...
% 	'NextPlot' , 'add' )
% axis vis3d
% patch( ...
% 	struct( ...
% 	'faces' , fv.tri.faces , ...
% 	'vertices' , fv.tri.vertices ) , ...
% 	'FaceAlpha' , 1.0 , ...
% 	'FaceColor' , [1,0.7,0.7] , ...
% 	'EdgeColor' , [0,0,0] , ...
% 	'EdgeAlpha' , 1 )
% for ff = 1 : size(fv.tri.faces,1)
% 	vv = fv.tri.faces(ff,:) ;
% 	x = fv.tri.vertices(vv,:) ;
% 	c = mean( x ) ;
% 	n = cross( x(2,:)-x(1,:) , x(3,:)-x(2,:) ) ;
% 	n = n / norm(n) ;
% 	n = n * 0.1*a ;
% 	plot3( ...
% 		c(1) + [0,n(1)] , ...
% 		c(2) + [0,n(2)] , ...
% 		c(3) + [0,n(3)] , ...
% 		'bo-' , ...
% 		'MarkerFaceColor' , 'b' , ...
% 		'LineWidth' , 3 )
% 	if( dot( n , c ) <= 0 )
% 		disp( ff )
% 	end
% end
% text( ...
% 	1.1 * fv.tri.vertices(:,1) , ...
% 	1.1 * fv.tri.vertices(:,2) , ...
% 	1.1 * fv.tri.vertices(:,3) , ...
% 	cellstr(num2str((1:size(fv.tri.vertices,1))')) , ...
% 	'FontSize' , 20 , ...
% 	'FontWeight' , 'bold' )

%% add thick walls

h = 0.10 * a ;
nv = size( fv.tri.vertices ,1) ;
nn = nan( [size(fv.tri.faces,1),1] ) ;

% triangular finite element with rectangular walls
fv.tri.vertices = [
	fv.tri.vertices
	fv.tri.vertices * (1-(h/a))
	] ;
faces = [ fv.tri.faces , fv.tri.faces+nv , nn ] ;
ind = [
	1,2,3,7
	6,5,4,7
	1,4,5,2
	1,3,6,4
	2,5,6,3
	] ;
fv.tri.faces = [] ;
for ii = 1 : size(faces,1)
	tmp = faces(ii,:) ;
	fv.tri.faces = [
		fv.tri.faces
		tmp(ind)
		] ;
end
fv.tri.elements = reshape( 1:size(fv.tri.faces,1) ,5,[])' ;

%% VIEW

% figure( ...
% 	'Color' , [1,1,1] )
% axes( ...
% 	'DataAspectRatio' , [1,1,1] , ...
% 	'Visible' , 'off' , ...
% 	'NextPlot' , 'add' )
% axis vis3d
% patch( ...
% 	struct( ...
% 	'faces' , fv.tri.faces , ...
% 	'vertices' , fv.tri.vertices ) , ...
% 	'FaceAlpha' , 1.0 , ...
% 	'FaceColor' , [1,0.7,0.7] , ...
% 	'EdgeColor' , [0,0,0] , ...
% 	'EdgeAlpha' , 1 )
% for ee = 1 : size( fv.tri.elements ,1)
% 	for ff = 1 : 2
% 		vv = fv.tri.faces( fv.tri.elements(ee,ff) ,1:3) ;
% 		x = fv.tri.vertices(vv,:) ;
% 		c = mean( x ) ;
% 		n = cross( x(2,:)-x(1,:) , x(3,:)-x(2,:) ) ;
% 		n = n / norm(n) ;
% 		n = n * 0.1*a ;
% 		plot3( ...
% 			c(1) + [0,n(1)] , ...
% 			c(2) + [0,n(2)] , ...
% 			c(3) + [0,n(3)] , ...
% 			'bo-' , ...
% 			'MarkerFaceColor' , 'b' , ...
% 			'LineWidth' , 3 )
% 		switch ff
% 			case 1
% 				if( dot( n , c ) <= 0 )
% 					disp( ff )
% 				end
% 			case 2
% 				if( dot( n , c ) >= 0 )
% 					disp( ff )
% 				end
% 		end
% 	end
% end
% text( ...
% 	1.1 * fv.tri.vertices(:,1) , ...
% 	1.1 * fv.tri.vertices(:,2) , ...
% 	1.1 * fv.tri.vertices(:,3) , ...
% 	cellstr(num2str((1:size(fv.tri.vertices,1))')) , ...
% 	'FontSize' , 20 , ...
% 	'FontWeight' , 'bold' )

%%

end

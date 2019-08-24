
function [ nodes , springs ] = make_geom_truncOct_3D( geom_size )
% 3D soccer ball

% geom_size = [3,3,3] ;
% tissue = number of alveoli in block [ #x , #y , #z ]

rotAng = [ 0 , pi/4 , 0 ] ; % radians, rotate geometry about x-axis, y-axis, z-axis
geomType = 'tissue' ;
B2 = constructGeometry( geomType , geom_size , rotAng ) ;
v = B2.vertices ;
edges = [
	B2.faces(:,1) , B2.faces(:,2)
	B2.faces(:,2) , B2.faces(:,3)
	B2.faces(:,3) , B2.faces(:,1)
	] ;
edges = sort( edges ,2) ;
edges = unique( edges ,'rows') ;
edges_length = sqrt(sum(power( v(edges(:,2),:) - v(edges(:,1),:) ,2),2)) ;
clear B2

%
num_points = size(v,1) ;
num_springs = size(edges,1) ;
num_dimensions = size(v,2) ;
num_stiffness_tension = 1 ;
num_stiffness_compression = 0 ;

%
nodes.position = v ;
nodes.fixed = false([num_points,1]) ;
nodes.force = zeros( size(nodes.position) ) ;
springs.nodes = edges ;
springs.stiffness_tension     = zeros([num_springs,num_stiffness_tension    ]) + 1.0 ;
springs.stiffness_compression = zeros([num_springs,num_stiffness_compression]) + 1.0 ;
springs.restlength = 0.25 * edges_length ;
springs.compression = false([num_springs,1]) ;

% boundary conditions
p_bound = 30*[
	2
	2
	1
	1
	3
	3
	] ;
dir_bound = {
	[-1, 0,0]
	[+1, 0,0]
	[ 0,-1,0]
	[ 0,+1,0]
	[ 0,0,-1]
	[ 0,0,+1]
	} ;
ind_bound_v = {
	find( abs( v(:,1) - min(v(:,1)) ) / range(v(:,1)) < 1e-2 )
	find( abs( v(:,1) - max(v(:,1)) ) / range(v(:,1)) < 1e-2 )
	find( abs( v(:,2) - min(v(:,2)) ) / range(v(:,2)) < 1e-2 )
	find( abs( v(:,2) - max(v(:,2)) ) / range(v(:,2)) < 1e-2 )
	find( abs( v(:,3) - min(v(:,3)) ) / range(v(:,3)) < 1e-2 )
	find( abs( v(:,3) - max(v(:,3)) ) / range(v(:,3)) < 1e-2 )
	} ;

for ii = 1 : numel(ind_bound_v)
	nodes.force( ind_bound_v{ii} ,:) = nodes.force( ind_bound_v{ii} ,:) + (dir_bound{ii}*p_bound(ii)/numel(ind_bound_v{ii})) ;
end

% figure( ...
% 	'Color' , [1,1,1] )
% style = { 'g-' , 'c-' , 'r:' , 'b:' } ;
% axes( ...'XLim' , [-5,10] , ...'YLim' , [-5,10] , ...'ZLim' , [-5,10] , ...
% 	'DataAspectRatio' , [1,1,1] , ...
% 	'NextPlot' , 'add' )
% for ii = 1 : num_springs
% 	plot3( ...
% 		nodes.position(springs.nodes(ii,:),1) , ...
% 		nodes.position(springs.nodes(ii,:),2) , ...
% 		nodes.position(springs.nodes(ii,:),3) , ...
% 		'k-' , ...
% 		'LineWidth' , 1 )
% end
% for ii = 1 : num_points
% 	plot3( ...
% 		nodes.position(ii,1) + [0,nodes.force(ii,1)] , ...
% 		nodes.position(ii,2) + [0,nodes.force(ii,2)] , ...
% 		nodes.position(ii,3) + [0,nodes.force(ii,3)] , ...
% 		'r-' , ...
% 		'LineWidth' , 1 )
% end

end

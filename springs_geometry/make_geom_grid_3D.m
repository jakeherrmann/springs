function [ nodes , springs ] = make_geom_grid_3D( geom_size )

%%

[ x , y , z ] = meshgrid( ...
	0 : geom_size(1) , ...
	0 : geom_size(2) , ...
	0 : geom_size(3) ) ;
v = [ x(:) , y(:) , z(:) ] ;
ind = v + 1 ;
neighborhood = [
	-1 ,  0 ,  0
	+1 ,  0 ,  0
	 0 , -1 ,  0
	 0 , +1 ,  0
	 0 ,  0 , -1
	 0 ,  0 , +1
	] ;
N = size(v,1) ;
M = size(neighborhood,1) ;
ind = bsxfun(@plus, ind , permute(neighborhood,[3,2,1]) ) ;
ind( (ind<1) ) = nan ;
for ii = 1 : 3
	for kk = 1 : M
		remove = ind(:,ii,kk) > size(x,ii) ;
		ind(remove,ii,kk) = nan ;
	end
end
ind = sub2ind( size(x) , ind(:,2,:) , ind(:,1,:) , ind(:,3,:) ) ;
ind = permute( ind ,[1,3,2]) ;
ind = cat(3, ind , repmat((1:N)',[1,M]) ) ;
ind = reshape( ind , [ N*M , 2 ] ) ;
ind( any(isnan(ind),2) ,:) = [] ;
ind = unique( sort(ind,2) ,'rows') ;

%%

num_springs = size(ind,1) ;
num_stiffness_tension = 1 ;
num_stiffness_compression = 0 ;

nodes.position = v ;
nodes.fixed = false([N,1]) ;
nodes.force = zeros( size(nodes.position) ) ;

springs.nodes = ind ;
springs.stiffness_tension     = zeros([num_springs,num_stiffness_tension    ]) + 1.0 ;
springs.stiffness_compression = zeros([num_springs,num_stiffness_compression]) + 1.0 ;
springs.rest_length = 0.5 + zeros([num_springs,1]) ;
springs.compression = false([num_springs,1]) ;

springs.stiffness_tension = normrnd( 1.0 , 0.5 , [num_springs,num_stiffness_tension] ) ;
springs.stiffness_tension( springs.stiffness_tension < 0.1 ) = 1.0 ;

springs.rest_length = springs.rest_length + max(0, normrnd(0.0,0.3,size(springs.rest_length)) ) ;

%%

% boundary conditions
p_bound = 10*[
	0
	0
	1
	1
	1
	1
	] ;
fix_bound = [
	true
	true
	false
	false
	false
	false
	] ;

p_bound = 10*[
	0
	0
	0
	0
	0
	0
	] ;
fix_bound = [
	true
	true
	true
	true
	true
	true
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
	nodes.fixed( ind_bound_v{ii} ) = nodes.fixed( ind_bound_v{ii} ) | fix_bound(ii) ;
end

%%

end
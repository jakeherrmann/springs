function [ nodes , springs ] = make_geom_hexagon_2D( geom_size )
% 2D HEXAGONAL

% geom_size = [ 10 , 7 ] ;

num_row = geom_size(1) + 4 ;
num_col = geom_size(2) + 3 ;

[ X , Y ] = meshgrid( ...
	0 : (num_row-1) , ...
	0 : (num_col-1) ) ;
Y = bsxfun(@plus, Y , 0.5*mod(1:num_row,2) ) ;
X = [ X(:) , Y(:) ] ;

[ v , c ] = voronoin( X ) ;
ind_remove = find( (v(:,2)<0) | (v(:,2)>(num_col-1)) ) ;
c( cellfun( @(cc)any(ismember(cc,ind_remove)) , c ) ) = [] ;
c = cellfun( @(cc)cc([1:end,1]) , c , 'UniformOutput' , false ) ;

c = cellfun( @(cc) [ cc(1:(end-1)) ; cc(2:end) ]' , c , 'UniformOutput' , false ) ;
c = cell2mat(c) ;
c = sort( c ,2) ;
c = unique( c ,'rows') ;

ind_remove = find( v(:,1) == 0.3750 ) ;
c( any(ismember(c,ind_remove),2) ,:) = [] ;

ind_remove = find( v(:,1) == (num_row-2)+0.6250 ) ;
c( any(ismember(c,ind_remove),2) ,:) = [] ;

ind_remove = find( ( v(c(:,1),2) == 0.5000 ) & ( v(c(:,2),2) == 0.5000 ) ) ;
c( ind_remove ,:) = [] ;

ind_remove = find( ( v(c(:,1),2) == (num_col-1) ) & ( v(c(:,2),2) == (num_col-1) ) ) ;
c( ind_remove ,:) = [] ;

p_bound = [
	1
	1
	1
	1
	] ;
dir_bound = {
	[-1, 0]
	[+1, 0]
	[ 0,-1]
	[ 0,+1]
	} ;
ind_bound_v = {
	find( ( v(:,1) == 0.6250 ) )
	find( ( v(:,1) == (num_row-2)+0.3750 ) )
	find( ( v(:,2) == 0.5000 ) )
	find( ( v(:,2) == (num_col-1) ) )
	} ;
ind_bound = any(ismember( c , cat(1,ind_bound_v{:}) ),2) ;
b = c( ind_bound ,:) ;
c( ind_bound ,:) = [] ;
for ii = 1 : size(b,1)
	b(ii,:) = b(ii, 1+ismember( b(ii,:) , cat(1,ind_bound_v{:}) ) ) ;
end
d = cell(size(ind_bound_v)) ;
for ii = 1 : numel(d)
	d{ii} = find( ismember( b(:,2) , ind_bound_v{ii} ) ) ;
end

v(:,1) = v(:,1) * (0.3125/0.375) ;

f = zeros( [ size(b,1) , size(v,2) ] ) ;
for ii = 1 : numel(d)
	n = v( b(d{ii},2) ,:) - v( b(d{ii},1) ,:) ;
	n = bsxfun(@rdivide, n , sqrt(sum(power( n ,2),2)) ) ;
	n_dot_dir = sum(bsxfun(@times, n , dir_bound{ii} ),2) ;
	f( d{ii} ,:) = n * ( sum(n_dot_dir) * p_bound(ii) ) / numel(d{ii}) ;
end

nodes_id = unique( c(:) ) ;
nodes.fixed    = false( size(nodes_id) ) ;
nodes.position = v( nodes_id ,:) ;
nodes.force    = zeros( size(nodes.position) ) ;
for ii = 1 : size(b,1)
	ind = find( nodes_id == b(ii,1) ) ;
	nodes.force(ind,:) = nodes.force(ind,:) + f(ii,:) ;
end
num_points = numel( nodes_id ) ;

num_stiffness_tension = 1 ;
num_stiffness_compression = 0 ;

num_springs = size(c,1) ;
[ ~ , springs.nodes ] = ismember( c , nodes_id ) ;
springs.stiffness_tension     = zeros([num_springs,num_stiffness_tension    ]) + 1.0 ;
springs.stiffness_compression = zeros([num_springs,num_stiffness_compression]) + 1.0 ;
springs.restlength  = zeros([num_springs,1]) + 0.25 ;
springs.compression = false([num_springs,1]) ;

% figure( ...
% 	'Color' , [1,1,1] )
% style = { 'g-' , 'c-' , 'r:' , 'b:' } ;
% axes( ...
% 	'XLim' , [-5,15] , ...
% 	'YLim' , [-5,15] , ...
% 	'DataAspectRatio' , [1,1,1] , ...
% 	'NextPlot' , 'add' )
% for ii = 1 : size(c,1)
% 	plot( ...
% 		v(c(ii,:),1) , ...
% 		v(c(ii,:),2) , ...
% 		'k-' , ...
% 		'LineWidth' , 2 )
% end
% for ii = 1 : numel(d)
% 	for jj = 1 : numel(d{ii})
% 		plot( ...
% 			v(b(d{ii}(jj),1),1) + [0,f(d{ii}(jj),1)] , ...
% 			v(b(d{ii}(jj),1),2) + [0,f(d{ii}(jj),2)] , ...
% 			style{ii} , ...
% 			'LineWidth' , 2 )
% 	end
% end

end

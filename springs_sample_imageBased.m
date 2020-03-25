%% DESCRIPTION

% this example shows a rectangular grid in 2D, with
% spring properties and fixed boundaries derived from
% an image. spring compression is not allowed, and all
% springs are initially pre-stressed (length greater
% than rest length).  Good solver algorithm/objective
% pairs for this problem are newton/energy and
% steepest/energy.  steepest/sumforce also seems to do
% okay.

%%

addpath( 'springs_matlab_solver' )
addpath( 'springs_matlab_visualizer' )

%% parameters for rectangular grid

% create a phantom image
img = phantom( 'Modified Shepp-Logan' , 64 ) ;
ind = (img<=0) & imfill(img>0,'holes') ;
img(ind) = 0.1 ;
[ num_row , num_col ] = size( img ) ;

% convert image value to desired spring constant [0,1] and rest length
img_k = zeros(size(img)) ;
img_k(img>0) = img(img>0)*2.5 ;
img_L0 = 0.4 + zeros(size(img)) ;

%% PREPROCESS

% row and column index of each node in the rectangular grid
[ col , row ] = meshgrid( 1:num_col , 1:num_row ) ;

% location of each node, each row of this array corresponds
% to one nodes's row index and column index
loc = cat(2, ...
	reshape( row ,[],1) , ...
	reshape( col ,[],1) ) ;

% adjacent nodes are identified from each node looking backwards
% to the row-1 and col-1 indexes, forming a 2D lattice with
% 4-point connectivity.
adj = [
	bsxfun(@plus, loc , [ 0,-1] )
	bsxfun(@plus, loc , [-1, 0] )
	] ;

% each spring connects two adjacent nodes.
% each row of this array correspond to one spring's start and
% end nodes.
spring_nodes = cat(2, ...
	repmat( (1:size(loc,1))' ,[2,1]) , ...
	adj(:,1) + num_row*(adj(:,2)-1) ) ;

% remove springs with nodes that point outside the boundary
remove = any( bsxfun(@lt,adj,[1,1]) | bsxfun(@gt,adj,[num_row,num_col]) ,2) ;
spring_nodes( remove ,:) = [] ;

% use image value to determine spring stiffnesses and rest lengths.
% each spring represents the average of the two voxels it connects.
spring_stiffness   = mean( img_k (spring_nodes) ,2) ;
spring_rest_length = mean( img_L0(spring_nodes) ,2) ;

% remove springs with zero stiffness.
remove = spring_stiffness <= 0 ;
spring_stiffness( remove ) = [] ;
spring_rest_length( remove ,:) = [] ;
spring_nodes( remove ,:) = [] ;

% remove very stiff springs, and replace by boundary/fixed nodes.
remove = spring_stiffness > 0.999 ;
boundary = {
	reshape( spring_nodes( remove ,:) ,[],1)
	} ;
spring_stiffness( remove ,:) = [] ;
spring_rest_length( remove ,:) = [] ;
spring_nodes( remove ,:) = [] ;

% % alternatively, could fix nodes around the image boundary
% boundary = {
% 	find(ismember(loc(:,1),1      ))
% 	find(ismember(loc(:,1),num_row))
% 	find(ismember(loc(:,2),1      ))
% 	find(ismember(loc(:,2),num_col))
% 	} ;
% boundary{3} = setdiff( boundary{3} , [boundary{1};boundary{2}] ) ;
% boundary{4} = setdiff( boundary{4} , [boundary{1};boundary{2}] ) ;

%% CREATE NODES/SPRINGS DATA STRUCTURES

% overall parameters
num_dimensions = 2 ;            % number of spatial dimensions
precision = 'double' ;          % floating point precision: 'single' or 'double'
num_stiffness_tension = 1 ;     % number of polynomial stiffness coefficients for tension
num_stiffness_compression = 0 ; % number of polynomial stiffness coefficients for compression

% size of network
num_points = size(loc,1) ;
num_springs = size(spring_nodes,1) ;

% these are the only mandatory properties of nodes data structure.
% the springs_solve() command will ignore any other properties.
% you may add any additional properties as desired.
nodes.position = zeros( [ num_points , num_dimensions ] ) ;
nodes.force    = zeros( [ num_points , num_dimensions ] ) ;
nodes.fixed    = false( [ num_points , num_dimensions ] ) ;

% assign values computed in previous section
nodes.position = [ col(:) , row(:) ] ;
for bb = 1 : numel(boundary)
	nodes.fixed( boundary{bb} ,:) = true ;
end

% these are the only mandatory properties of springs data structure.
% the springs_solve() command will ignore any other properties.
% you may add any additional properties as desired.
springs.nodes                 = zeros( [ num_springs , 2                         ] ) ;
springs.rest_length           = zeros( [ num_springs , 1                         ] ) ;
springs.stiffness_tension     = zeros( [ num_springs , num_stiffness_tension     ] ) ;
springs.stiffness_compression = zeros( [ num_springs , num_stiffness_compression ] ) ;
springs.compression           = zeros( [ num_springs , 1                         ] ) ;

% assign values computed in previous section
springs.nodes = spring_nodes ;
springs.rest_length(:) = spring_rest_length ;
springs.stiffness_tension(:,1) = spring_stiffness ;

%% SOLVE EQUILIBRIUM (NO FORCES)

% solver parameters not specified are treated as default.
options = springs_default_options() ;
options.algorithm = 'steepest' ;
options.objective = 'energy' ;
options.include_force_fixed_nodes = true ;
options.num_iter_print = 1000 ;
options.tolerance_change_objective = 1e-24 ;
options.tolerance_sum_net_force = 1e-24 ;

% run solver
[ nodes_eq , ~ ] = springs_solve( nodes , springs , options ) ;

%% DISPLACE THE FIXED BOUNDARY OUTWARDS AND RESOLVE EQUILIBRIUM

% select a few random nodes (not on boundary) and apply random forces
xy0 = 0.5*size(img) ;
ind = all( nodes_eq.fixed ,2) ;
xy = nodes_eq.position(ind,:) ;
xy = bsxfun(@plus, xy , -xy0 ) ;
xy = xy * 1.2 ;
xy = bsxfun(@plus, xy , +xy0 ) ;
nodes_eq_stretch = nodes_eq ;
nodes_eq_stretch.position(ind,:) = xy ;

% run solver
[ nodes_eq_stretch , ~ ] = springs_solve( nodes_eq_stretch , springs , options ) ;

%% DISPLAY

h_before = display_2D( nodes_eq         , springs , 'stiffness_tension' , [0,1] , parula(10) , true ) ;
h_after  = display_2D( nodes_eq_stretch , springs , 'stiffness_tension' , [0,1] , parula(10) , true ) ;
linkaxes( [h_before.ax,h_after.ax] , 'xy' ) ;

%%


%% DESCRIPTION

% this example shows a rectangular grid in 2D, with all fixed
% boundaries, and external forces applied at random locations.
% spring compression is allowed, and all springs are initially
% at their rest lengths. best solver algorithm/objective seems
% to be anneal/maxforce for this scenario.

%%

addpath( 'springs_matlab_solver' )
addpath( 'springs_matlab_visualizer' )

%% parameters for rectangular grid

num_row = 20 ;
num_col = 15 ;
k_mean = 1.0 ;
k_stdv = 0.0 ;
L0_mean = 1.0 ;
L0_stdv = 0.2 ;

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

% randomly sample a normal distribution of spring stiffnesses
spring_stiffness = normrnd( k_mean , k_stdv , [size(spring_nodes,1),1] ) ;

% randomly sample a normal distribution of spring rest lengths
spring_rest_length = normrnd( L0_mean , L0_stdv , [size(spring_nodes,1),1] ) ;

% remove springs with zero stiffness.
remove = spring_stiffness <= 0 ;
spring_stiffness( remove ) = [] ;
spring_rest_length( remove ,:) = [] ;
spring_nodes( remove ,:) = [] ;

% fix nodes around the boundary
boundary = {
	find(ismember(loc(:,1),1      ))
	find(ismember(loc(:,1),num_row))
	find(ismember(loc(:,2),1      ))
	find(ismember(loc(:,2),num_col))
	} ;
boundary{3} = setdiff( boundary{3} , [boundary{1};boundary{2}] ) ;
boundary{4} = setdiff( boundary{4} , [boundary{1};boundary{2}] ) ;

%% CREATE NODES/SPRINGS DATA STRUCTURES

% overall parameters
num_dimensions = 2 ;            % number of spatial dimensions
precision = 'double' ;          % floating point precision: 'single' or 'double'
num_stiffness_tension = 1 ;     % number of polynomial stiffness coefficients for tension
num_stiffness_compression = 1 ; % number of polynomial stiffness coefficients for compression

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
nodes.fixed( boundary{1} ,:) = true ;
nodes.fixed( boundary{2} ,:) = true ;
nodes.fixed( boundary{3} ,:) = true ;
nodes.fixed( boundary{4} ,:) = true ;

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
springs.stiffness_compression(:,1) = spring_stiffness ;
springs.compression(:) = true ;

%% SOLVE EQUILIBRIUM (NO FORCES)

% solver parameters not specified are treated as default.
options = springs_default_options() ;
options.algorithm = 'anneal' ;
options.objective = 'maxforce' ;
options.include_force_fixed_nodes = false ;
options.num_iter_print = 10000 ;
options.tolerance_change_objective = 1e-24 ;
options.tolerance_sum_net_force = 1e-24 ;

% run solver
[ nodes_eq , ~ ] = springs_solve( nodes , springs , options ) ;

%% ADD FORCES AND RESOLVE EQUILIBRIUM

% select a few random nodes (not on boundary) and apply random forces
num_rand_force = 10 ;
magnitude = 0.2 ;
free_node_ind = find( ~any( nodes_eq.fixed ,2) ) ;
rand_node_ind = free_node_ind( randperm( numel(free_node_ind) ,num_rand_force) ) ;
rand_force_ang = rand( [num_rand_force,1] ) *(2.0*pi) ;
rand_force_mag = zeros( size(rand_force_ang) ) + magnitude ;
nodes_eq.force( rand_node_ind ,:) = bsxfun(@times, rand_force_mag , [cos(rand_force_ang),sin(rand_force_ang)] ) ;

% run solver
[ nodes_eq_force , ~ ] = springs_solve( nodes_eq , springs , options ) ;

%% DISPLAY

h_before = display_2D( nodes_eq       , springs , 'stiffness_tension' , [0,2] , parula(10) , true ) ;
h_after  = display_2D( nodes_eq_force , springs , 'stiffness_tension' , [0,2] , parula(10) , true ) ;
linkaxes( [h_before.ax,h_after.ax] , 'xy' ) ;

%%


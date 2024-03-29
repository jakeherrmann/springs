%% DESCRIPTION

% this example shows a cubic grid in 3D, with all fixed
% boundaries, and external forces applied at random locations.
% spring compression is allowed, and all springs are initially
% at their rest lengths. best solver algorithm/objective seems
% to be anneal/maxforce for this scenario.

%%

addpath( 'springs_matlab_solver' )
addpath( 'springs_matlab_visualizer' )

%% parameters for rectangular grid

num_row = 7 ;
num_col = 8 ;
num_lay = 9 ;
k_mean = 1.0 ;
k_stdv = 0.1 ;
L0_mean = 1.0 ;
L0_stdv = 0.02 ;

%% PREPROCESS

% row, column, layer index of each node in the rectangular grid
[ col , row , lay ] = meshgrid( 1:num_col , 1:num_row , 1:num_lay ) ;

% location of each node, each row of this array corresponds
% to one nodes's row index and column index
loc = cat(2, ...
	reshape( row ,[],1) , ...
	reshape( col ,[],1) , ...
	reshape( lay ,[],1) ) ;

% adjacent nodes are identified from each node looking backwards
% to the row-1, col-1, lay-1 indexes, forming a 3D lattice with
% 6-point connectivity.
adj = [
	bsxfun(@plus, loc , [ 0,-1, 0] )
	bsxfun(@plus, loc , [-1, 0, 0] )
	bsxfun(@plus, loc , [ 0, 0,-1] )
	] ;

% each spring connects two adjacent nodes.
% each row of this array correspond to one spring's start and
% end nodes.
spring_nodes = cat(2, ...
	repmat( (1:size(loc,1))' ,[3,1]) , ...
	adj(:,1) + num_row*(adj(:,2)-1) + num_row*num_col*(adj(:,3)-1) ) ;

% remove springs with nodes that point outside the boundary
remove = any( bsxfun(@lt,adj,[1,1,1]) | bsxfun(@gt,adj,[num_row,num_col,num_lay]) ,2) ;
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
	find(ismember(loc(:,3),1      ))
	find(ismember(loc(:,3),num_lay))
	} ;
boundary{3} = setdiff( boundary{3} , [boundary{1};boundary{2}] ) ;
boundary{4} = setdiff( boundary{4} , [boundary{1};boundary{2}] ) ;
boundary{5} = setdiff( boundary{5} , [boundary{1};boundary{2};boundary{3};boundary{4}] ) ;
boundary{6} = setdiff( boundary{6} , [boundary{1};boundary{2};boundary{3};boundary{4}] ) ;

%% CREATE NODES/SPRINGS DATA STRUCTURES

% overall parameters
num_dimensions = 3 ;            % number of spatial dimensions
precision = 'double' ;          % floating point precision: 'single' or 'double'

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
nodes.position = [ col(:) , row(:) , lay(:) ] ;
nodes.fixed( boundary{1} ,:) = true ;
nodes.fixed( boundary{2} ,:) = true ;
nodes.fixed( boundary{3} ,:) = true ;
nodes.fixed( boundary{4} ,:) = true ;
nodes.fixed( boundary{5} ,:) = true ;
nodes.fixed( boundary{6} ,:) = true ;

% these are the only mandatory properties of springs data structure.
% the springs_solve() command will ignore any other properties.
% you may add any additional properties as desired.
springs.nodes                               = zeros( [ num_springs , 2 ] ) ;
springs.rest_length                         = zeros( [ num_springs , 1 ] ) ;
springs.force_length_type_tension           = zeros( [ num_springs , 1 ] ) ;
springs.force_length_type_compression       = zeros( [ num_springs , 1 ] ) ;
springs.force_length_parameters_tension     =  cell( [ num_springs , 1 ] ) ;
springs.force_length_parameters_compression =  cell( [ num_springs , 1 ] ) ;

% assign values computed in previous section
springs.nodes = spring_nodes ;
springs.rest_length(:) = spring_rest_length ;
springs.force_length_type_tension(:) = 1 ; % none=0, polynomial=1, exponential=2, powerlaw=3
springs.force_length_parameters_tension = mat2cell( spring_stiffness , ones([num_springs,1]) , size(spring_stiffness,2) ) ;
springs.force_length_type_compression(:) = 1 ; % none=0, polynomial=1, exponential=2, powerlaw=3
springs.force_length_parameters_compression = mat2cell( spring_stiffness , ones([num_springs,1]) , size(spring_stiffness,2) ) ;

%% SOLVE EQUILIBRIUM (NO FORCES)

% solver parameters not specified are treated as default.
options = springs_default_options() ;
options.algorithm = 'steepest' ;
options.objective = 'maxforce' ;
options.include_force_fixed_nodes = false ;
options.num_iter_print = 10000 ;
options.tolerance_change_objective = 1e-12 ;
options.tolerance_sum_net_force = 1e-24 ;
options.use_numerical_hessian = true ;

% run solver
[ nodes_eq , ~ ] = springs_solve( nodes , springs , options ) ;

%% ADD FORCES AND RESOLVE EQUILIBRIUM

% select a few random nodes (not on boundary) and apply random forces
num_rand_force = 20 ;
magnitude = 1.0 ;
free_node_ind = find( ~any( nodes_eq.fixed ,2) ) ;
rand_node_ind = free_node_ind( randperm( numel(free_node_ind) ,num_rand_force) ) ;
rand_force_ang1 = rand( [num_rand_force,1] ) *(2.0*pi) ;
rand_force_ang2 = rand( [num_rand_force,1] ) * pi ;
rand_force_mag = zeros( size(rand_force_ang1) ) + magnitude ;
nodes_eq.force( rand_node_ind ,:) = bsxfun(@times, ...
	rand_force_mag , ...
	[sin(rand_force_ang2).*cos(rand_force_ang1),sin(rand_force_ang2).*sin(rand_force_ang1),cos(rand_force_ang2)] ) ;

% run solver
[ nodes_eq_force , ~ ] = springs_solve( nodes_eq , springs , options ) ;

%% DISPLAY

h_before = display_3D( nodes_eq       , springs , springs_tension(nodes_eq      ,springs) , [0,0.2] , parula(10) , true ) ;
h_after  = display_3D( nodes_eq_force , springs , springs_tension(nodes_eq_force,springs) , [0,0.2] , parula(10) , true ) ;
linkprop( [h_before.ax,h_after.ax] , {'XLim','YLim','ZLim','View'} ) ;

%%

%% DESCRIPTION

% this example shows a rectangular grid in 2D, with partially fixed
% boundaries, and a gravitational field of external forces applied
% everywhere. spring compression is not allowed, and all springs
% are initially prestressed (length greater than rest length).
% best solver algorithm/objective seems to be steepest/sumforce.
% anneal/sumforce is okay but takes much longer.

%%

addpath( 'springs_matlab_solver' )
addpath( 'springs_matlab_visualizer' )

%% parameters for rectangular grid

num_row = 20 ;
num_col = 18 ;
k_mean = 1.0 ;
k_stdv = 0.2 ;
L0_mean = 0.2 ;
L0_stdv = 0.0 ;

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
nodes.fixed( boundary{2} ,:) = true ; % fix top row vertically
nodes.fixed( boundary{1} ,:) = true ; % bottom row fixed
nodes.fixed( boundary{3} ,1) = true ; % fix left column horizontally
nodes.fixed( boundary{4} ,1) = true ; % fix right column horizontally

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
options.objective = 'sumforce' ;
options.include_force_fixed_nodes = true ;
options.num_iter_print = 1000 ;
options.tolerance_change_objective = 1e-24 ;
options.tolerance_sum_net_force = 1e-24 ;

% run solver
[ nodes_eq , ~ ] = springs_solve( nodes , springs , options ) ;

%% ADD FORCES AND RESOLVE EQUILIBRIUM

% use distance from the "top" in a gravitational field
gravity_angle = -0.35*pi ;
gravity_direction = [ cos(gravity_angle) , sin(gravity_angle) ] ;
height = nodes.position * gravity_direction' ;
gravity_magnitude = 0.2 * ( max(height) - height ) / range(height) ;
nodes_eq.force = bsxfun(@times, gravity_magnitude , gravity_direction ) ;

% run solver
[ nodes_eq_force , ~ ] = springs_solve( nodes_eq , springs , options ) ;

%% DISPLAY

h_before = display_2D( nodes_eq       , springs , 'stiffness_tension' , [0,2] , parula(10) , true ) ;
h_after  = display_2D( nodes_eq_force , springs , 'stiffness_tension' , [0,2] , parula(10) , true ) ;
linkaxes( [h_before.ax,h_after.ax] , 'xy' ) ;

%%


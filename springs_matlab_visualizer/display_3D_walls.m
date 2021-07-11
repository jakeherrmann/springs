function h = display_3D( nodes , springs , walls , varargin )
% display_3D( nodes , springs ) plots all springs as black edges.
% 
% h = display_3D( nodes , springs ) returns handles to the new figure and axes.
%
% display_3D( nodes , springs , color_variable , color_range ) plots all springs with color
% determined according to color_variable linearly within min/max of color_range.  Color_variable
% may be either a string corresponding to a field of springs data structure, or it may be an array
% of size NS x 1 where NS is the number of springs. Default color_range is [0,1].  Default variable
% is none.
%
% display_3D( nodes , springs , color_variable , color_range , color_map ) plots all springs
% with the number of colors NC, and RGB values of each color, specified as an NC x 3 array
% color_map. Default is jet(10)
%
% display_3D( nodes , springs , color_variable , color_range , color_map , show_forces ) also plots
% the external applied force at each node with black arrows, scaled by the largest force magnitude,
% if show_forces is true.  Default is false.
%
% Example:
% display_3D( nodes , springs , 'stiffness_tension' , [0,1] , jet(11) , true )
%


% default parameters
color_variable = '' ;
color_range = [0,1] ;
color_map = jet(10) ;
show_forces = false ;

% user-defined parameters
if numel(varargin) > 0
	color_variable = varargin{1} ;
end
if numel(varargin) > 1
	color_range = varargin{2} ;
	color_range = sort( color_range ) ;
end
if numel(varargin) > 2
	color_map = varargin{3} ;
end
if numel(varargin) > 3
	show_forces = varargin{4} ;
end

% group springs according to color_variable
% or make all springs black by default
if isempty( color_variable )
	color_map = [0,0,0] ;
	color_ind = ones( [ size(springs.nodes,1) , 1 ] ) ;
elseif ischar( color_variable )
	if ~isfield( springs , color_variable )
		fprintf( 'ERROR: %s is not a property of springs\n' , color_variable )
		color_map = [0,0,0] ;
		color_ind = ones( [ size(springs.nodes,1) , 1 ] ) ;
	else
		color_val = springs.(color_variable)(:,1) ;
		color_ind = (color_val-min(color_range))/range(color_range) ;
		color_ind = min(1, max(0, color_ind )) ;
		color_ind = 1 + round( (size(color_map,1)-1) * color_ind ) ;
	end
else
	if ~isequal( size(color_variable,1) , size(springs.nodes,1) )
		color_map = [0,0,0] ;
		color_ind = ones( [ size(springs.nodes,1) , 1 ] ) ;
	else
		color_val = color_variable ;
		color_ind = (color_val-min(color_range))/range(color_range) ;
		color_ind = min(1, max(0, color_ind )) ;
		color_ind = 1 + round( (size(color_map,1)-1) * color_ind ) ;
	end
end

%
h.fig = figure( ...
	'Colormap' , color_map , ...
	'Color' , [1,1,1] ) ;
h.ax = axes( ...'XLim' , [ -1 10 ] , ...'YLim' , [ -1 10 ] , ...'ZLim' , [ -1 10 ] , ...
	'CLim' , color_range , ...
	'View' , [ -25 , 25 ] , ...
	'DataAspectRatio' , [1,1,1] , ...
	'Clipping' , 'off' , ...
	'TickDir' , 'out' , ...
	'NextPlot' , 'add' ) ;

fv.vertices = nodes.position ;
fv.faces = cellfun( @(w) unique( reshape( springs.nodes( w ,:)' ,[],1) ,'stable')' , walls ,'UniformOutput',false) ;
flip = cellfun( @(w) ~ismember( springs.nodes(w(1),1) , springs.nodes(w(end),:) ) , walls ) ;
max_wall_size = max(cellfun(@numel,fv.faces)) ;
for ii = 1 : numel(fv.faces)
	if flip(ii)
		fv.faces{ii}([1,2]) = fv.faces{ii}([2,1]) ;
	end
	fv.faces{ii}( (end+1) : (max_wall_size+1) ) = nan ;
end
fv.faces = cell2mat( fv.faces ) ;
patch( fv , ...
	'AmbientStrength'  , 0.2 , ...
	'DiffuseStrength'  , 0.6 , ...
	'SpecularStrength' , 0.0 , ...
	'FaceColor' , [0,0,0]+0.7 , ...
	'FaceAlpha' , 0.3 , ...
	'EdgeColor' , 'none' )

for ii = 1 : size(color_map,1)
	ind = find( color_ind == ii ) ;
	if ~isempty( ind )
		xx = [ nodes.position( springs.nodes(ind,1) ,1) , nodes.position( springs.nodes(ind,2) ,1) ] ;
		yy = [ nodes.position( springs.nodes(ind,1) ,2) , nodes.position( springs.nodes(ind,2) ,2) ] ;
		zz = [ nodes.position( springs.nodes(ind,1) ,3) , nodes.position( springs.nodes(ind,2) ,3) ] ;
		xx(:,3) = nan ;
		yy(:,3) = nan ;
		zz(:,3) = nan ;
		plot3( ...
			reshape( xx' ,[],1) , ...
			reshape( yy' ,[],1) , ...
			reshape( zz' ,[],1) , ...
			'-' , ...
			'Color' , color_map(ii,:) , ...
			'LineWidth' , 1 ) ;
	end
end

axis tight
axis_size = max( [ range(xlim()) , range(ylim()) , range(zlim()) ] ) ;

if show_forces
	force_mag_sq = sum(power( nodes.force ,2),2) ;
	ind = (force_mag_sq>0) & (force_mag_sq>(1e-3*max(force_mag_sq))) ;
	force_scaled = nodes.force(ind,:) / sqrt(max(force_mag_sq)) ;
	
	force_scaled = force_scaled ;
	
	xx = bsxfun(@plus, nodes.position(ind,1) , bsxfun(@times, force_scaled(:,1) , [0,1] ) ) ;
	yy = bsxfun(@plus, nodes.position(ind,2) , bsxfun(@times, force_scaled(:,2) , [0,1] ) ) ;
	zz = bsxfun(@plus, nodes.position(ind,3) , bsxfun(@times, force_scaled(:,3) , [0,1] ) ) ;
	xx(:,end+1) = nan ;
	yy(:,end+1) = nan ;
	zz(:,end+1) = nan ;
	
	plot3( ...
		reshape( xx' ,[],1) , ...
		reshape( yy' ,[],1) , ...
		reshape( zz' ,[],1) , ...
		'k-' , ...
		'LineWidth' , 3 )
end

axis tight

if nargout > 0
	varargout{1} = h ;
end

end
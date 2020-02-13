function varargout = display_2D( nodes , springs , varargin )
% display_2D( nodes , springs ) plots all springs as black edges.
% 
% h = display_2D( nodes , springs ) returns handles to the new figure and axes.
%
% display_2D( nodes , springs , color_variable , color_range ) plots all springs with color
% determined according to color_variable linearly within min/max of color_range.  Color_variable
% may be either a string corresponding to a field of springs data structure, or it may be an array
% of size NS x 1 where NS is the number of springs.
%
% display_2D( nodes , springs , color_variable , color_range , color_map ) plots all springs
% with the number of colors NC, and RGB values of each color, specified as an NC x 3 array
% color_map.
%
% Example:
% display_2D( nodes , springs , 'stiffness_tension' , [0,1] , jet(11) )
%

% default parameters
color_variable = '' ;
color_range = [0,1] ;
color_map = jet(10) ;

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

h.fig = figure( ...
	'Color' , [1,1,1] ) ;
h.ax = axes( ...
	'DataAspectRatio' , [1,1,1] , ...
	'TickDir' , 'out' , ...
	'NextPlot' , 'add' ) ;

for ii = 1 : size(color_map,1)
	ind = find( color_ind == ii ) ;
	if ~isempty( ind )
		xx = [ nodes.position( springs.nodes(ind,1) ,1) , nodes.position( springs.nodes(ind,2) ,1) ] ;
		yy = [ nodes.position( springs.nodes(ind,1) ,2) , nodes.position( springs.nodes(ind,2) ,2) ] ;
		xx(:,3) = nan ;
		yy(:,3) = nan ;
		plot( ...
			reshape( xx' ,[],1) , ...
			reshape( yy' ,[],1) , ...
			'-' , ...
			'Color' , color_map(ii,:) , ...
			'LineWidth' , 2 ) ;
	end
end

axis tight

if nargout > 0
	varargout{1} = h ;
end

end
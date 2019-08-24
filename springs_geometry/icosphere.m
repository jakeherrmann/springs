function [vv,ff] = icosphere(varargin)
%ICOSPHERE Generate icosphere.
% Create a unit geodesic sphere created by subdividing a regular
% icosahedron with normalised vertices.
%
%   [V,F] = ICOSPHERE(N) generates to matrices containing vertex and face
%   data so that patch('Faces',F,'Vertices',V) produces a unit icosphere
%   with N subdivisions.
%
%   FV = ICOSPHERE(N) generates an FV structure for using with patch.
%
%   ICOSPHERE(N) and just ICOSPHERE display the icosphere as a patch on the
%   current axes and does not return anything.
%
%   ICOSPHERE uses N = 3.
%
%   ICOSPHERE(AX,...) plots into AX instead of GCA.
%
%   See also SPHERE.
%
%   Based on C# code by Andres Kahler
%   http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

% Parse possible axes input
if nargin > 2
	error('Too many input variables, must be 0, 1 or 2.');
end
[cax,args,nargs] = axescheck(varargin{:});

n = 3; % default number of sub-divisions
if nargs > 0, n = args{1}; end % override based on input

% generate regular unit icosahedron (20 faced polyhedron)
[v,f] = icosahedron(); % size(v) = [12,3]; size(f) = [20,3];

% ALLOCATE VERTICES AND FACES
nv = 12 ;
v( 2*12*power(4,n) ,:) = [0,0,0] ;

% recursively subdivide triangle faces
for gen = 1:n
	
	% for each triangle
	indf = 1 : (20*power(4,gen-1)) ;
	indv = permute( f(indf,:) , [2,3,1] ) ;
	
	% calculate mid points (add new points to v)
	newv = 0.5*( v(circshift(indv,-1),:) + v(indv,:) ) ;
	newv = bsxfun(@rdivide, newv , sqrt(sum(power(newv,2),2)) ) ;
	noldv = nv ;
	nnewv = size(newv,1) ;
	
	%
	v( noldv + (1:nnewv) ,:) = newv ;
	nv = noldv + nnewv ;
	
	% generate new subdivision triangles
	newf = cat(2, ...
		indv , ...
		bsxfun(@plus, [1,3;2,1;3,2] , permute( noldv + 3*(0:(size(indv,3)-1)) ,[3,1,2]) ) ) ;
	newf = [
		reshape( permute( newf , [2,1,3] ) ,3,[],1).'
		bsxfun(@plus, noldv + 3*(0:(size(indv,3)-1)).' , [1,2,3] )
		] ;
	
	% replace triangle with subdivision
	f = newf ;
end
v = v(1:nv,:) ;

% remove duplicate vertices
[v,b,ix] = unique(v,'rows') ;
% reassign faces to trimmed vertex list and remove any duplicate faces
f = unique(ix(f),'rows');

switch(nargout)
	case {0,1} % return fv structure for patch
		vv = struct('Vertices',v,'Faces',f,...
			'VertexNormals',v,'FaceVertexCData',v(:,3));
	case 2 % return vertices and faces
		vv = v; ff = f;
	otherwise
		error('Too many output variables, must be 0, 1 or 2.');
end

end

function [v,f] = icosahedron()
%ICOSAHEDRON creates unit regular icosahedron
%   Returns 12 vertex and 20 face values.
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; % v1
	1, t, 0; % v2
	-1,-t, 0; % v3
	1,-t, 0; % v4
	0,-1, t; % v5
	0, 1, t; % v6
	0,-1,-t; % v7
	0, 1,-t; % v8
	t, 0,-1; % v9
	t, 0, 1; % v10
	-t, 0,-1; % v11
	-t, 0, 1];% v12
% normalise vertices to unit size
v = v./norm(v(1,:));

% create faces
f = [ 1,12, 6; % f1
	1, 6, 2; % f2
	1, 2, 8; % f3
	1, 8,11; % f4
	1,11,12; % f5
	2, 6,10; % f6
	6,12, 5; % f7
	12,11, 3; % f8
	11, 8, 7; % f9
	8, 2, 9; % f10
	4,10, 5; % f11
	4, 5, 3; % f12
	4, 3, 7; % f13
	4, 7, 9; % f14
	4, 9,10; % f15
	5,10, 6; % f16
	3, 5,12; % f17
	7, 3,11; % f18
	9, 7, 8; % f19
	10, 9, 2];% f20
end

function B = constructGeometry( geomType , geomSize , refinement , fe_type , rigidRotation )

%%

switch geomType
	case 'comsolCoarse'
		%% TEST: use exported comsol geometry
		B = readComsolMesh( fullfile( 'comsol_geometry' , 'coarse.mphtxt' ) ) ;
		
	case 'comsolFine'
		%% TEST: use exported comsol geometry
		B = readComsolMesh( fullfile( 'comsol_geometry' , 'fine.mphtxt' ) ) ;
	
	case 'tissue'
		%% construct a 3D tessellated block of alveoli
		
		nx = geomSize(1) ;
		ny = geomSize(2) ;
		nz = geomSize(3) ;
		
		latticeType = 1 ;
		
		switch latticeType
			case 1
				% grow 1 row
				xx = (1:(nx-1))' ;
				oo = zeros([nx-1,1]) ;
				conn = [ xx+1 , oo+4 , xx , oo+1 , oo+4 ] ;
				
				% add on rows
				xx = (1:nx)' ;
				oo = zeros([nx,1]) ;
				for yy = 2 : ny
					switch mod(yy,2)
						case 0; ff = [5,1] ;
						case 1; ff = [8,4] ;
					end
					conn = [
						conn
						xx+nx*(yy-1) , oo+ff(1) , xx+nx*(yy-2) , oo+ff(2) , oo+6
						] ;
				end
				
				% add on layers
				for zz = 2 : nz
					for yy = 1 : ny
						switch mod(yy,2)
							case 0
								switch mod(zz,2)
									case 0; ff = [3,7] ;
									case 1; ff = [6,2] ;
								end
							case 1
								switch mod(zz,2)
									case 0; ff = [6,2] ;
									case 1; ff = [3,7] ;
								end
						end
						conn = [
							conn
							xx+nx*(yy-1)+nx*ny*(zz-1) , oo+ff(1) , xx+nx*(yy-1)+nx*ny*(zz-2) , oo+ff(2) , oo+6
							] ;
					end
				end
				
				rotAng = [
					-0.5*( pi - acos(-1/3) )
					0.0
					0.0
					] ;
				
			case 2
				% grow 1 row
				xx = (1:(nx-1))' ;
				oo = zeros([nx-1,1]) ;
				conn = [ xx+1 , oo+6 , xx , oo+2 , oo+6 ] ;
				
				% add on rows
				xx = (1:nx)' ;
				oo = zeros([nx,1]) ;
				for yy = 2 : ny
					switch mod(yy,2)
						case 0; ff = [5,1] ;
						case 1; ff = [8,4] ;
					end
					conn = [
						conn
						xx+nx*(yy-1) , oo+ff(1) , xx+nx*(yy-2) , oo+ff(2) , oo+6
						] ;
				end
				
				% add on layers
				for zz = 2 : nz
					for yy = 1 : ny
						switch mod(zz,3)
							case 0; ff = [3,7] ; ffn = 6 ;
							case 2; ff = [3,7] ; ffn = 6 ;
							case 1; ff = [1,4] ; ffn = 4 ;
						end
						conn = [
							conn
							xx+nx*(yy-1)+nx*ny*(zz-1) , oo+ff(1) , xx+nx*(yy-1)+nx*ny*(zz-2) , oo+ff(2) , oo+ffn
							] ;
					end
				end
				
				rotAng = [
					0.0
					0.0
					0.0
					] ;
		end
		
		% rotate to align with coordinate axes
		rotMatrix = rotmat( rotAng ) ;
		numAlv = nx * ny * nz ;
		for ii = numAlv : -1 : 1
			a(ii) = alvGeom() ;
			a(ii).sqr.vertices = ( rotMatrix * a(ii).sqr.vertices.' ).' ;
			a(ii).hex.vertices = ( rotMatrix * a(ii).hex.vertices.' ).' ;
			a(ii).tri.vertices = ( rotMatrix * a(ii).tri.vertices.' ).' ;
		end
		
		% centroid positions
		c = zeros( [ numel(a) , 3 ] ) ;
		for ii = 1 : size(conn,1)
			jj  = conn(ii,1) ;
			jjf = conn(ii,2) ;
			kk  = conn(ii,3) ;
			kkf = conn(ii,4) ;
			ffn = conn(ii,5) ;
			switch ffn
				case 4; ss = 'sqr' ;
				case 6; ss = 'hex' ;
			end
			dc = diff( [
				a(jj).(ss).vertices( a(jj).(ss).faces(jjf,1) ,:)
				a(kk).(ss).vertices( a(kk).(ss).faces(kkf,1) ,:)
				] ,1,1) ;
			c(jj,:) = c(kk,:) + dc ;
		end
		for ii = 1 : numel(a)
			for ss = {'sqr','hex','tri'}
				a(ii).(ss{1}).vertices = bsxfun(@plus, a(ii).(ss{1}).vertices , c(ii,:) ) ;
			end
		end
		
		%
		N = size( a(1).sqr.vertices ,1) ;
		for ss = {'sqr','hex'}
			b = [a.(ss{1})] ;
			bb.(ss{1}) = cell2mat( arrayfun( @(i) b(i).faces + (i-1)*N ,1:numel(b) ,'UniformOutput' ,false)' ) ;
		end
		B.nodes = cat(1,b.vertices) ;
		B.springs = [
			reshape( permute( cat(3, bb.sqr , circshift( bb.sqr ,[0,-1]) ) ,[3,2,1]) ,2,[])'
			reshape( permute( cat(3, bb.hex , circshift( bb.hex ,[0,-1]) ) ,[3,2,1]) ,2,[])'
			] ;
		B.walls = [
			mat2cell( reshape( (1:numel(bb.sqr)) + 0             ,size(bb.sqr,2),[])' ,ones(1,size(bb.sqr,1)),size(bb.sqr,2))
			mat2cell( reshape( (1:numel(bb.hex)) + numel(bb.sqr) ,size(bb.hex,2),[])' ,ones(1,size(bb.hex,1)),size(bb.hex,2))
			] ;
		[ B.nodes , ~ , ind ] = uniquetol( B.nodes ,1e-6,'ByRows',true) ;
		B.springs = ind(B.springs) ;
		[ B.springs , ~ , ind ] = unique( sort( B.springs ,2) ,'rows') ;
		B.walls = cellfun( @(w) ind(w)' , B.walls ,'UniformOutput',false) ;
		
	case 'cube'
		%% TEST: make a cube out of 2 triangular prisms
		
		vertices = [
			-1 , -1 , -1
			+1 , -1 , -1
			-1 , +1 , -1
			-1 , -1 , +1
			+1 , +1 , -1
			+1 , -1 , +1
			-1 , +1 , +1
			+1 , +1 , +1
			] ;
		vertices = bsxfun(@times, vertices , geomSize(:)' ) ;
		alignFaces = [
			1 , 2 , 3 , 0
			4 , 5 , 6 , 0
			1 , 6 , 5 , 2
			1 , 3 , 4 , 6
			2 , 5 , 4 , 3
			] ;
		
		faces = [] ;
		facesTmp = zeros(size(alignFaces)) ;
		for jj = 1 : 2
			switch jj
				case 1, indTmp = [7,4,6,2,1,3] ;
				case 2, indTmp = [6,8,7,3,5,2] ;
			end
			for ii = 1 : numel(facesTmp)
				if alignFaces(ii)==0
					facesTmp(ii) = nan ;
				else
					facesTmp(ii) = indTmp(alignFaces(ii)) ;
				end
			end
			faces = [
				faces
				facesTmp
				] ;
		end
		elements = [
			1:5
			6:10
			] ;
		B = struct( ...
			'elements' , elements , ...
			'faces'    , faces    , ...
			'vertices' , vertices ) ;
		
	case 'prism'
		%% TEST: a single triangular prism
		
		vertices = [
			0 , 1 , +1
			0 , 0 , +1
			1 , 0 , +1
			1 , 0 , -1
			0 , 0 , -1
			0 , 1 , -1
			] ;
		vertices = bsxfun(@times, vertices , geomSize(:)' ) ;
		alignFaces = [
			1 , 2 , 3 , 0
			4 , 5 , 6 , 0
			1 , 6 , 5 , 2
			1 , 3 , 4 , 6
			2 , 5 , 4 , 3
			] ;
		
		faces = alignFaces ;
		faces( faces==0 ) = nan ;
		elements = [
			1:5
			] ;
		B = struct( ...
			'elements' , elements , ...
			'faces'    , faces    , ...
			'vertices' , vertices ) ;
		
	case 'ball'
		%% TEST: make a ball from a tesselated icosahedron
		
		% project tesselated icosahedron onto sphere
		fv = icosphere(2) ;
		fv = struct( ...
			'vertices' , fv.Vertices , ...
			'faces' , fv.Faces ) ;
		
		% resize radius
		a = 1.0 ;
		fv.vertices = fv.vertices * a ;
		
		% add thick walls
		h = 0.10 * a ;
		nv = size( fv.vertices ,1) ;
		nn = nan( [size(fv.faces,1),1] ) ;
		
		fv.vertices = [
			fv.vertices
			fv.vertices * (1-(h/a))
			] ;
		faces = [ fv.faces , fv.faces+nv , nn ] ;
		ind = [
			1,2,3,7
			6,5,4,7
			1,4,5,2
			1,3,6,4
			2,5,6,3
			] ;
		B.vertices = fv.vertices ;
		B.faces = [] ;
		for ii = 1 : size(faces,1)
			tmp = faces(ii,:) ;
			B.faces = [
				B.faces
				tmp(ind)
				] ;
		end
		B.elements = reshape( 1:size(B.faces,1) ,5,[])' ;
		
end

%%

% rigidly rotate geometry
R = rotmat(rigidRotation) ;
B.nodes = B.nodes * R' ;

%%

end
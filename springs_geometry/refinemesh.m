function B = refinemesh( B , edgeType , numDiv )

if numDiv<=1
	return
end

weight = (1:(numDiv-1)) / numDiv ;
weight = permute( weight ,[3,1,2]) ;
weight = cat(1, 1-weight , weight ) ;

basis = {} ;
vertices = B.vertices ;
faces = {} ;
indNodes = [1,6,11,2,7,12] ;

switch edgeType
	case 'tri'
		if numDiv > 2
			error( 'ERROR: refinemesh: Splitting triangular edges more than twice not yet defined' )
		end
		for ee = 1 : size(B.elements,1)
			oldfaces = B.faces( B.elements(ee,:) ,:) ;
			newvertices = nan([6*(numDiv-1),3]) ;
			indNewVertices = (1:size(newvertices,1)) + size(vertices,1) ;
			indNewVertices = reshape( indNewVertices ,[numDiv-1,6]) ;
			for jj = 1 : 6
				tt  = ceil(jj/3) ;
				ii1 = 1 + mod(jj-1,3) ;
				ii2 = 1 + mod(jj  ,3) ;
				if( tt==2 )
					ii1 = 1 + mod(-ii1,3) ;
					ii2 = 1 + mod(-ii2,3) ;
				end
				vv = [ oldfaces(tt,ii1) , oldfaces(tt,ii2) ] ;
				xx = sum( bsxfun(@times, vertices(vv,:) , weight ) ,1) ;
				xx = permute( xx , [3,2,1] ) ;
				newvertices( (1:(numDiv-1))+(jj-1)*(numDiv-1) ,:) = xx ;
			end
			vertices = [
				vertices
				newvertices
				] ;
			% the following only works for numDiv=2
			indNewVertices = reshape( indNewVertices ,[3,2]) ;
			indNewVertices = [
				[ oldfaces(1,[1,2,3])' , oldfaces(2,[3,2,1])' ]
				indNewVertices
				] ;
			indTri = [
				1,4,6
				5,6,4
				4,2,5
				6,5,3
				] ;
			[a,b] = ismember( oldfaces , oldfaces(indNodes) ) ;
			newfaces = repmat( {oldfaces} ,[size(indTri,1),1]) ;
			for kk = 1 : size(indTri,1)
				c = indNewVertices( indTri(kk,:) ,:) ;
				c = cat(2, ...
					c([1,2,3],1) , ...
					c([3,2,1],2) ) ;
				newfaces{kk}(a) = c(b(a)) ;
			end
			faces = [
				faces
				newfaces
				] ;
			basis = [
				basis
				repmat( B.elementBasis(ee) ,[size(indTri,1),1])
				] ;
		end
		
	case 'rect'
		for ee = 1 : size(B.elements,1)
			oldfaces = B.faces( B.elements(ee,:) ,:) ;
			newvertices = nan([3*(numDiv-1),3]) ;
			indNewVertices = (1:size(newvertices,1)) + size(vertices,1) ;
			indNewVertices = reshape( indNewVertices ,[numDiv-1,3]) ;
			ind1 = [1,2,3] ;
			ind2 = [3,2,1] ;
			for jj = 1:3
				ii1 = ind1(jj) ;
				ii2 = ind2(jj) ;
				vv = [ oldfaces(1,ii1) , oldfaces(2,ii2) ] ;
				xx = sum( bsxfun(@times, vertices(vv,:) , weight ) ,1) ;
				xx = permute( xx , [3,2,1] ) ;
				newvertices( (1:(numDiv-1))+(jj-1)*(numDiv-1) ,:) = xx ;
			end
			vertices = [
				vertices
				newvertices
				] ;
			indNewVertices = [
				oldfaces(1,ind1)
				indNewVertices
				oldfaces(2,ind2)
				] ;
			[a,b] = ismember( oldfaces , oldfaces(indNodes) ) ;
			newfaces = repmat( {oldfaces} ,[numDiv,1]) ;
			for kk = 1 : numDiv
				c = cat(2, ...
					indNewVertices(kk  ,ind1) , ...
					indNewVertices(kk+1,ind2) ) ;
				newfaces{kk}(a) = c(b(a)) ;
			end
			faces = [
				faces
				newfaces
				] ;
			basis = [
				basis
				repmat( B.elementBasis(ee) ,[numDiv,1])
				] ;
		end
end

faces = cat(1,faces{:}) ;

% recompose the geometry structure
B.vertices = vertices ;
B.faces = faces ;
B.elements = reshape( 1:size(faces,1) ,5,[])' ;
B.elementBasis = basis ;

B = removeDuplicateVertices( B ) ;

end
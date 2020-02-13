function force = springs_tension( nodes , springs )
% force = springs_tension( nodes , springs ) returns the force in each spring, positive for
% tension and negative for compression, if compression is allowed.

delta_position = nodes.position( springs.nodes(:,2) ,:) - nodes.position( springs.nodes(:,1) ,:) ;
length = sqrt(sum(power( delta_position ,2),2)) ;
delta_length = length - springs.rest_length ;
force = zeros( [size(springs.nodes,1),1] ) ;

% springs in tension
NST = size( springs.stiffness_tension ,2) ;
if NST > 0
	ind = find( delta_length > 0.0 ) ;
	force(ind,:) = sum( bsxfun(@power,delta_length(ind),1:NST) .* springs.stiffness_tension(ind,:) ,2) ;
end

% springs in compression
NSC = size( springs.stiffness_compression ,2) ;
if NSC > 0
	ind = find( ( delta_length < 0.0 ) & springs.compression ) ;
	force(ind,:) = -sum( bsxfun(@power,delta_length(ind),1:NST) .* springs.stiffness_compression(ind,:) ,2) ;
end

end
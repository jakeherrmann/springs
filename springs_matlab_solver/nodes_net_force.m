function net_force = nodes_net_force( nodes , springs )
% force = springs_tension( nodes , springs ) returns the net force at each node.

springs_force = springs_tension( nodes , springs ) ;
springs_vector = nodes.position( springs.nodes(:,2) ,:) - nodes.position( springs.nodes(:,1) ,:) ;
springs_vector = bsxfun(@rdivide, springs_vector , sqrt(sum(power( springs_vector ,2),2)) ) ;
springs_force = bsxfun(@times, springs_force , springs_vector ) ;

% start with external forces, and add contribution from each spring
net_force = nodes.force ;
for dd = 1 : size(nodes.position,2)
	net_force(:,dd) = net_force(:,dd) + accumarray( ...
		[  springs.nodes(:,1)  ;  springs.nodes(:,2)  ] , ...
		[ +springs_force(:,dd) ; -springs_force(:,dd) ] , ...
		[size(net_force,1),1] , ...
		@sum , ...
		0 ) ;
end

end


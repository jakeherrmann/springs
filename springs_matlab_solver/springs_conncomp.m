function springs_label = springs_conncomp( nodes , springs )
	
	num_springs = size(springs.nodes,1) ;
	springs_label = zeros( [num_springs,1] ) ;
	unlabeled = (1:num_springs)' ;
	nodes_fixed = find(all(nodes.fixed,2)) ;
	
	next_label = 1 ;
	while ~isempty( unlabeled )
		next_unlabeled = 1 ;
		while ~isempty( next_unlabeled )
			frontier = unlabeled( next_unlabeled ) ;
			springs_label(frontier) = next_label ;
			unlabeled(next_unlabeled) = [] ;
			nodes_curr = springs.nodes(frontier,:) ;
			nodes_curr = setdiff( nodes_curr , nodes_fixed ) ;
			next_unlabeled = find( any( ismember( springs.nodes(unlabeled,:) , nodes_curr ) ,2) ) ;
		end
		next_label = next_label + 1 ;
	end
	
end

function force = springs_tension( nodes , springs )
% force = springs_tension( nodes , springs ) returns the force in each spring, positive for
% tension and negative for compression, if compression is allowed.

num_springs = size(springs.nodes,1) ;

length = springs_length( nodes , springs ) ;
delta_length = length - springs.rest_length ;
force = zeros( [num_springs,1] ) ;
for ss = 1 : num_springs
	if delta_length(ss) == 0.0
		% spring at rest
		continue ;
	elseif delta_length(ss) > 0.0
		% spring in tension
		DL = delta_length(ss) ;
		force_length_type = springs.force_length_type_tension(ss) ;
		force_length_parameters = springs.force_length_parameters_tension{ss}(:)' ;
		force_sign = +1 ;
	elseif delta_length(ss) < 0.0
		% spring in compression
		DL = abs(delta_length(ss)) ;
		force_length_type = springs.force_length_type_compression(ss) ;
		force_length_parameters = springs.force_length_parameters_compression{ss}(:)' ;
		force_sign = -1 ;
	end
	switch force_length_type
		case 0 %none
			continue ;
		case 1 %polynomial
			force_magnitude = sum( bsxfun(@power,DL,1:numel(force_length_parameters)) .* force_length_parameters ,2) ;
		case 2 %exponential
			force_magnitude = sum( exp( bsxfun(@times,DL,force_length_parameters(2:2:end)) ) .* force_length_parameters(1:2:end) ,2) ;
		case 3 %powerlaw
			force_magnitude = sum( bsxfun(@power,DL,force_length_parameters(2:2:end)) .* force_length_parameters(1:2:end) ,2) ;
	end
	force(ss) = force_magnitude * force_sign ;
end

end
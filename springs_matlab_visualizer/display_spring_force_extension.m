function h = display_spring_force_extension( nodes , springs , strain )
% display_spring_force_extension( nodes , springs , max_strain )

num_springs = size(springs.nodes,1) ;

Lr = springs.rest_length ;
Lr0 = Lr ; % initial rest length

strain = strain(:)' ;
stress = nan( num_springs , numel(strain) ) ;
for ee = 1 : numel(strain)
	length = (1.0+strain(ee)) * Lr0 ;
	delta_length = length - Lr ;
	%
	failed = (length./Lr0) > springs.failure_threshold ;
	%
	yielding = (length./Lr) > springs.plastic_threshold ;
	modifier = ( ((length./Lr)./springs.plastic_threshold).*springs.plastic_modifier + (1-springs.plastic_modifier) ) ;
	Lr(yielding) = Lr(yielding) .* modifier(yielding) ;
	%
	force = zeros( [num_springs,1] ) ;
	for ss = 1 : num_springs
		if failed(ss)
			force(ss) = 0.0 ;
			continue
		end
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
	stress(:,ee) = force ;
end

strain = repmat( [ strain , strain(end) , nan ] ,[num_springs,1]) ;
stress = [ stress , zeros(size(stress,1),1) , nan(size(stress,1),1) ] ;

%%

scaling = 'log' ;

h.fig = figure( ...
	'Color' , [1,1,1] ) ;
h.ax = axes( ...
	'XScale' , scaling , ...
	'YScale' , scaling , ...
	'NextPlot' , 'add' ) ;
h.patch = patch( ...
	100 * reshape( strain' ,[],1) , ...
	max(1e-4, reshape( stress' ,[],1) ) , ...
	'k' , ...
	'FaceColor' , 'none' , ...
	'EdgeColor' , [0,0,0] , ...
	'EdgeAlpha' , 0.1 ) ;
xlabel( 'Strain (%)' )
ylabel( 'Stress' )

%%

end
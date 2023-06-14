function length = springs_length( nodes , springs )
% length = springs_length( nodes , springs ) returns the length of each spring

delta_position = nodes.position( springs.nodes(:,2) ,:) - nodes.position( springs.nodes(:,1) ,:) ;
length = sqrt(sum(power( delta_position ,2),2)) ;

end
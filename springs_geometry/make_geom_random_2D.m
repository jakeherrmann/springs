function [ nodes , springs ] = make_geom_random_2D()

%%

N=20;
seed=12;
ratio=1/3;
pCross = 1 ;

[p,adj,collagen,elastin]=generateIsotropicNetwork(N,seed,pCross,ratio);
adj=logical(adj);

edgeNodes=sortNodes(p);
internalNodes=setdiff(1:size(p,1), ...
    [edgeNodes{1}' edgeNodes{2}' edgeNodes{3}' edgeNodes{4}' ]);

%initialize isotropic appearence
[Aeq,beq]=genConstraint(length(adj),edgeNodes);
d=squareform(pdist(p));
r=rand(length(adj/2)).*d.*triu(adj);
r=r+r';
k=r*0;
k(adj)=1./r(adj);

%% set spring constants and equilibrium lengths
n=size(collagen,1);
d=squareform(pdist(p));

r=.95*d.*adj;
randR=1+betarnd(5,5,n,n);
randR=triu(randR)+triu(randR)';
r(collagen)=randR(collagen).*d(collagen);
r(elastin)=d(elastin);
% r(collagen)=d(collagen);


k=(1./d).*adj;
k(isnan(k))=0;
k(collagen)=d(collagen).^(-1);
k(elastin)=(d(elastin)).^(-1);
% scale=mean(k(elastin))/min(k(collagen));
% k(collagen)=k(collagen)*scale;
% bndry=setdiff(1:length(adj)^2, [collagen; elastin]);
% k(bndry)=k(bndry)/scale;

%colBefore{q}=r(collagen)./d(collagen); 

%%


num_points = size(p,1) ;
nodes.position = p ;
nodes.force    = zeros( size(p) ) ;
nodes.fixed    = false( size(p,1) ) ;
nodes.fixed( cat(1,edgeNodes{:}) ) = true ;

num_stiffness_tension = 1 ;
num_stiffness_compression = 0 ;

[node_start,node_end] = find( adj ) ;
springs.nodes = unique( sort( [ node_start , node_end ] ,2) ,'rows') ;
num_springs = size(springs.nodes,1) ;
springs.stiffness_tension     = zeros([num_springs,num_stiffness_tension    ]) + 1.0 ;
springs.stiffness_compression = zeros([num_springs,num_stiffness_compression]) + 1.0 ;
springs.restlength  = zeros([num_springs,1]) + 0.025 ;
springs.compression = false([num_springs,1]) ;

%%

end


function [nodes,N,collagen, elastin]=generateIsotropicNetwork(N,seed,pCross,ratio)

if seed>0
    rng(seed);
end
% pCross=1;
% N=1000;
% ratio=1/3;
% rng(12)
%% INITIALIZE MATRICES
%create empty matrices 
%(xBeg(n,1), yBeg(n,1)) connects to (xEnd(n,1), yEnd(n,1))
xBeg=zeros(N,1);
xEnd=xBeg;
yBeg=xBeg;
yEnd=xBeg;

randPoint=rand(N,2);

%% GENERATE LINES
    theta=rand(N,1)*2*pi;
    m=tan(theta);
    b=randPoint(:,2)-m.*randPoint(:,1);
%get boundary intersections
    x0_bound=b;
    x1_bound=m+b; 
    y0_bound=-b./m; 
    y1_bound=(1-b)./m; 
%find intersections that are in the unit square
    x0_idx=find(0<x0_bound & x0_bound<1);
    x1_idx=find(0<x1_bound & x1_bound<1);
    y0_idx=find(0<y0_bound & y0_bound<1);
    y1_idx=find(0<y1_bound & y1_bound<1);
% get boundary nodes
    x0_nodes=[0*x0_idx x0_bound(x0_idx)];
    x1_nodes=[0*x1_idx+1 x1_bound(x1_idx)];
    y0_nodes=[y0_bound(y0_idx) 0*y0_idx];
    y1_nodes=[y1_bound(y1_idx) 0*y1_idx+1];
    corner_nodes=[0 0; 0 1; 1 1; 1 0];
%combine matrices
    EdgeNodes=[x0_nodes; x1_nodes; y0_nodes; y1_nodes; corner_nodes];
    line_idx=[x0_idx; x1_idx; y0_idx; y1_idx];
%add internal lines
    for i=1:N
        line_coord=find(line_idx==i);
        xBeg(i)=EdgeNodes(line_coord(1),1);
        xEnd(i)=EdgeNodes(line_coord(2),1);
        yBeg(i)=EdgeNodes(line_coord(1),2);
        yEnd(i)=EdgeNodes(line_coord(2),2);
    end
%% GET INTERSECTIONS
XY=[xBeg yBeg xEnd yEnd];
XY=[corner_nodes [corner_nodes(2:4,:); corner_nodes(1,:)]; XY];
out = lineSegmentIntersect(XY,XY);

%% CROSSLINKING
adj=out.intAdjacencyMatrix;
inter=adj(5:end,5:end);
interNum=sum(inter(:))/2;

I=round(pCross*interNum);
nodeSet=randperm(interNum,I);
interIdx=find(triu(inter));
interNew=0*inter;
interNew(interIdx(nodeSet) )=1;
interNew=interNew+interNew';
adj(5:end,5:end)=interNew;

%% CREATE NODE ADJACENCY MATRIX
upTri=triu(adj);

matX=triu(out.intMatrixX,1);
matY=triu(out.intMatrixY,1);

matX(isnan(matX))=0;
matY(isnan(matY))=0;

nodes=[matX(upTri) matY(upTri)];
col=zeros(size(nodes,2),1);
row=col;
fiber=col;
n=0;

for i=1:N+4
    line=find(adj(i,:)==1);
    xLine=out.intMatrixX(i,line);
    yLine=out.intMatrixY(i,line);
    
    lineNodes=[xLine' yLine'];
    sortedMat=sortrows(lineNodes,[1 2]);
    [~, idx]=ismembertol(sortedMat,nodes,'ByRows',true);
    m=2*(length(idx)-1);
%     col(n+1:n+m)=[idx; circshift(idx,1)];
%     row(n+1:n+m)=[circshift(idx,1); idx];
    col(n+1:n+m)=[idx(1:end-1); idx(2:end)];
    row(n+1:n+m)=[idx(2:end); idx(1:end-1)];
    fiber(n+1:n+m)=i*ones(m,1);
    n=n+m;
end

[pairs,ac,~] = unique([col row],'rows');

adjFiber=sparse(pairs(:,1),pairs(:,2),fiber(ac));
A=logical(adjFiber);

%% SPLIT INTO FIBER TYPES
elasLines=randperm(N,round(ratio*N));
colLines=1:N;
colLines=setdiff(colLines, elasLines);
collagen=cell(length(colLines),1);
elastin=cell(length(elasLines),1);
for i=1:length(colLines)
    collagen{i}=find(adjFiber==4+colLines(i));
end
for i=1:length(elasLines)
    elastin{i}=find(adjFiber==4+elasLines(i));
end
% figure(6);
% G=graph(adjFiber);
% plot(G, 'XData',nodes(:,1),'YData',nodes(:,2));

function [p,adj,collagen, elastin]=generateIsotropicNetwork(N,seed,pCross,ratio)

% N=8;
% seed=12;
% pCross=1;
% ratio=1/3;
%% generateIsotropicNetwork: generate the random network
%
%% INPUTS:
%
% N == number of lines to initally generate
% seed == random seed for using repeated networks
% pCross == the proportion of crosslinking
%
%% OUTPUTS:
%
% xBeg == first x coordinate for each line
% xEnd == second x coordinate for each line
% yBeg == first y coordinate for each line
% yEnd == second y coordinate for each line
% nodes == [node number; x coordinate; y coordinate; a line index that ends at node;
%           a second line that ends there; a line index that begins at node;
%           a second line that begins there;  ]
% lengths == lengths of each line
%
%% LOGIC
%
% On the unit square, randomly seed N points
% For each point generate a random slope
% Draw a line through a point with the respective slope until it intersects
% the boundary
% Calculate all intersections and select those that occur within the unit
% square
% Bisect both lines at the point of intersection for each line
% Sort which lines intersect an edge and order them
% Add springs in between those points on the boundary
% A negative number (of the degree) is in place of would be lines for these nodes
%
if seed>0
    rng(seed);
end
% N=30;
% pCross=1;
%% INITIALIZE MATRICES
%create empty matrices 
%(xBeg(1,n), yBeg(1,n)) connects to (xEnd(1,n), yEnd(1,n)) creating a
%line
xBeg=zeros(1,N);
xEnd=zeros(1,N);
yBeg=zeros(1,N);
yEnd=zeros(1,N);
%list of nodes (1,:) are the x or y coordinates and (2:3,:) are the two
%lines at the intersection
xIntersection=zeros(3,1);
yIntersection=zeros(3,1);
%lengths and angles store just that which is indexed the same as the x and
%y matrices
randPoint=zeros(2,N);
% initially chose the fibers
elas=randperm(N,round(ratio*N));
col=1:N;
col=setdiff(col, elas);

%% GENERATE LINES
for n=1:N
    randPoint(1:2,n)=[rand, rand];
    for m=1:2
    used=false;
    theta=rand*pi;
        %check if it intersects x=0
    if (tan(theta)*(0-randPoint(1,n))+randPoint(2,n)<1) && (0<tan(theta)*(0-randPoint(1,n))+randPoint(2,n))
        x1=0;
        y1=tan(theta)*(0-randPoint(1,n))+randPoint(2,n);
        used=true;
    end
        %intersects x=1
    if (tan(theta)*(1-randPoint(1,n))+randPoint(2,n)<1) && (0<tan(theta)*(1-randPoint(1,n))+randPoint(2,n))
        if used==true
            x2=1;
            y2=tan(theta)*(1-randPoint(1,n))+randPoint(2,n);
        else
            x1=1;
            y1=tan(theta)*(1-randPoint(1,n))+randPoint(2,n);
            used=true;
        end
    end
    %intersects y=0
    if ((0-randPoint(2,n))/tan(theta)+randPoint(1,n)<1) && (0<(0-randPoint(2,n))/tan(theta)+randPoint(1,n)) 
        if used==true
            x2=(0-randPoint(2,n))/tan(theta)+randPoint(1,n);
            y2=0;
        else
            x1=(0-randPoint(2,n))/tan(theta)+randPoint(1,n);
            y1=0;
        end
    end
    %intersect y=1
    if ((1-randPoint(2,n))/tan(theta)+randPoint(1,n)<1) && (0<(1-randPoint(2,n))/tan(theta)+randPoint(1,n)) 
        x2=(1-randPoint(2,n))/tan(theta)+randPoint(1,n);
        y2=1;
    end
    %add new lines
    xBeg(1,n)=x1;    
    xEnd(1,n)=x2; 
    yBeg(1,n)=y1;    
    yEnd(1,n)=y2; 
    end
end

%% GET INTERSECTIONS
c=1;
for j=1:N %picks the first line
    x1Line=[xBeg(1,j) xEnd(1,j)];
    y1Line=[yBeg(1,j) yEnd(1,j)];
    for k=j+1:N %picks all lines coming after the first one
        x2Line=[xBeg(1,k) xEnd(1,k)];
        y2Line=[yBeg(1,k) yEnd(1,k)];
        %fit a line to the segment
        p1=polyfit(x1Line, y1Line,1);
        p2=polyfit(x2Line, y2Line,1);
        %get (x,y) of intersection 
        xIntersect = fzero(@(x) polyval(p1-p2,x),3); %set functions equal by letting their difference be zero and solve for x
        yIntersect = polyval(p1,xIntersect);
        %check if it is in the boundaries. If so, add it
        if ((0<=xIntersect) && (xIntersect<=1)) && ((0<=yIntersect) && (yIntersect<=1))
            xIntersection(1:3,c)=[xIntersect, j, k];
            yIntersection(1:3,c)=[yIntersect, j, k];
            c=c+1;
        end
    end
end
%% SPLIT LINES
numInter=size(xIntersection,2);
I=round(pCross*numInter);
nodesSet=randperm(numInter,I);
nodes=zeros(7,I);
index=1;

for m=nodesSet
%move coordinates and line numbers
    nodes(1,index)=m;
    nodes(2,index)=xIntersection(1,m);
    nodes(3,index)=yIntersection(1,m);
    %these lines end at the node
    nodes(4:5,index)=xIntersection(2:3,m);
    index=index+1;
end
for p=1:size(nodes,2)
    %new line number (add even)
    pE=(N-1)+(2*p);
    %new line number (add odd)
    pO=(N-1)+(2*p+1);
    %these lines begin at the original end point and end at the node
    nodes(6,p)=pE;
    nodes(7,p)=pO;
    
%move end points of the original lines so they are the beginning points of
%the new lines
    %line 3
    xEnd(pE)= xEnd(nodes(4,p));
    yEnd(pE)= yEnd(nodes(4,p));
    %line 4
    xEnd(pO)= xEnd(nodes(5,p));
    yEnd(pO)= yEnd(nodes(5,p));
%change original end points
%line 1
    xEnd(nodes(4,p))=nodes(2,p);
    yEnd(nodes(4,p))=nodes(3,p);
%line 2
    xEnd(nodes(5,p))=nodes(2,p);
    yEnd(nodes(5,p))=nodes(3,p);
%add beginning to new lines
%line 3
    xBeg(pE)=nodes(2,p);
    yBeg(pE)=nodes(3,p);
%line 4
    xBeg(pO)=nodes(2,p);
    yBeg(pO)=nodes(3,p);
% this monstrosity figures out what the others lines the two new lines
% intersect so this can be repeated
    for u=1:size(nodes,2)
        if u~=p %ignores the choosen intersection
            for v=4:5
                for w=4:5
                    if nodes(v,u) == nodes(w,p) %checks whether the lines in the intersection intersect any other lines
                        if xEnd(nodes(w,p))-xBeg(nodes(w,p))>0 %checks left to right orientation
                            if nodes(2,u) > nodes(2,p) %checks if the other intersection is on the right side of the line
                                if w==4
                                    nodes(v,u)= pE; %changes to the new line
                                else
                                    nodes(v,u)= pO; %changes to the new line
                                end
                            end
                        elseif xEnd(nodes(w,p))-xBeg(nodes(w,p))<0 %checks right to left orientation
                            if nodes(2,u) < nodes(2,p) %checks if the other intersection is on the left side of the line
                                if w==4
                                    nodes(v,u)= pE; %changes to the new line
                                else
                                    nodes(v,u)= pO; %changes to the new line
                                end
                            end
                        end
                    end
                end
            end
        end
    end
% keep new fibers their original type
line1=nodes(4,p);
line2=nodes(5,p);
line3=nodes(6,p);
line4=nodes(7,p);
%only need single comparison
if ismember(line1,col)==1
    col=[col line3];
else
    elas=[elas line3];
end
if ismember(line2,col)==1
    col=[col line4];
else
    elas=[elas line4];
end
    
end

% sort(col)

%% SORT
[rLines, lLines, tLines, bLines]=sortLines(xBeg, xEnd, yBeg, yEnd);
% put in ascending order
rLinesSort=sortrows(rLines);
lLinesSort=sortrows(lLines);
tLinesSort=sortrows(tLines);
bLinesSort=sortrows(bLines);

L=size(xBeg,2)+1;
M=size(nodes,2)+1;
M_0=max(max(nodes(4:7,:)))+1;
% create right side springs (x=1)
for i=1:(size(rLinesSort,1)+1)
    if i==1
        xBeg(1,L)=1;
        xEnd(1,L)=1;
        yBeg(1,L)=0;
        yEnd(1,L)=rLinesSort(i,1);
        nodes(1,M)=M;
        nodes(2,M)=1;
        nodes(3,M)=0;
        nodes(4,M)=-2;%change this one
        nodes(5,M)=-2;
        nodes(6,M)=L; 
        nodes(7,M)=-2;        
        N2=M; %mark corner node (bottom right)
        L=L+1;
        M=M+1;
    else
        xBeg(1,L)=1;
        xEnd(1,L)=1;
        yBeg(1,L)=yEnd(1,L-1);
        if i==size(rLinesSort,1)+1
            yEnd(1,L)=1;
        else
            yEnd(1,L)=rLinesSort(i,1);
        end
        nodes(1,M)=M;
        nodes(2,M)=1;
        nodes(3,M)=yBeg(1,L);
        nodes(4,M)=L-1;
        nodes(6,M)=L;
        if rLinesSort(i-1,3)==1
            nodes(7,M)=rLinesSort(i-1,2);
            nodes(5,M)=-3; %flag as degree 3 node
        else
            nodes(7,M)=-3; %flag as degree 3 node
            nodes(5,M)=rLinesSort(i-1,2);
        end
        E4=M;
        M=M+1;
        L=L+1;
    end
end
% create left side springs (x=1)
for i=1:(size(lLinesSort,1)+1)
    if i==1
        xBeg(1,L)=0;
        xEnd(1,L)=0;
        yBeg(1,L)=0;
        yEnd(1,L)=lLinesSort(i,1);
        nodes(1,M)=M;
        nodes(2,M)=0;
        nodes(3,M)=0;
        nodes(4,M)=-2;
        nodes(5,M)=-2;
        nodes(6,M)=L; 
        nodes(7,M)=-2; %change this one
        N1=M; % mark corner node (bottom left)
        L=L+1;
        M=M+1;
    else
        xBeg(1,L)=0;
        xEnd(1,L)=0;
        yBeg(1,L)=yEnd(1,L-1);
        if i==size(lLinesSort,1)+1
            yEnd(1,L)=1;
        else
            yEnd(1,L)=lLinesSort(i,1);
        end
        nodes(1,M)=M;
        nodes(2,M)=0;
        nodes(3,M)=yBeg(1,L);
        nodes(4,M)=L-1;
        nodes(6,M)=L;
        if lLinesSort(i-1,3)==1
            nodes(7,M)=lLinesSort(i-1,2);
            %flag as degree 3 node
            nodes(5,M)=-3;
        else
            %flag as degree 3 node
            nodes(7,M)=-3;
            nodes(5,M)=lLinesSort(i-1,2);
        end
        M=M+1;
        L=L+1;
    end
end
% create top springs
for i=1:(size(tLinesSort,1)+1)
    if i==1
        xBeg(1,L)=0;
        xEnd(1,L)=tLinesSort(i,1);
        yBeg(1,L)=1;
        yEnd(1,L)=1;
        nodes(1,M)=M;
        nodes(2,M)=0;
        nodes(3,M)=1;
        nodes(4,M)=nodes(6,M-1); %this finishes the top left corner node
        nodes(5,M)=-2;
        nodes(6,M)=L; 
        nodes(7,M)=-2;
        L=L+1;
        M=M+1;
    else
        yBeg(1,L)=1;
        yEnd(1,L)=1;
        xBeg(1,L)=xEnd(1,L-1);
        if i==size(tLinesSort,1)+1
            xEnd(1,L)=1;
        else
            xEnd(1,L)=tLinesSort(i,1);
        end
        nodes(1,M)=M;
        nodes(2,M)=xBeg(1,L);
        nodes(3,M)=1;
        nodes(4,M)=L-1;
        nodes(6,M)=L;
        if tLinesSort(i-1,3)==1
            nodes(7,M)=tLinesSort(i-1,2);
            %flag as degree 3 node
            nodes(5,M)=-3;
        else
            %flag as degree 3 node
            nodes(7,M)=-3;
            nodes(5,M)=tLinesSort(i-1,2);
        end
        M=M+1;
        L=L+1;
    end
end
%define top right corner node
nodes(1,M)=M;
nodes(2,M)=1;
nodes(3,M)=1;
nodes(4,M)=nodes(6,E4);
nodes(5,M)=L-1;
nodes(6,M)=-2;
nodes(7,M)=-2;
M=M+1;
%create bottom springs
for i=1:(size(bLinesSort,1)+1)
    if i==1
        xBeg(1,L)=0;
        xEnd(1,L)=bLinesSort(i,1);
        yBeg(1,L)=0;
        yEnd(1,L)=0;
        nodes(7,N1)=L;
        L=L+1;
    else
        yBeg(1,L)=0;
        yEnd(1,L)=0;
        xBeg(1,L)=xEnd(1,L-1);
        if i==size(bLinesSort,1)+1
            xEnd(1,L)=1;
        else
            xEnd(1,L)=bLinesSort(i,1);
        end
        nodes(1,M)=M;
        nodes(2,M)=xBeg(1,L);
        nodes(3,M)=0;
        nodes(4,M)=L-1;
        nodes(6,M)=L;
        if bLinesSort(i-1,3)==1
            nodes(7,M)=bLinesSort(i-1,2);
            %flag as degree 3 node
            nodes(5,M)=-3;
        else
            %flag as degree 3 node
            nodes(7,M)=-3;
            nodes(5,M)=bLinesSort(i-1,2);
        end
        M=M+1;
        L=L+1;
    end
end
nodes(4,N2)=L-1;

% define all edges as collagen

%col=[col M_0:max(max(nodes(6,:)))];

%% CONVERT TO ADJACENCY
nodes=nodes';
collagen=zeros(2*length(col),1);
indCol=1;
elastin=zeros(2*length(elas),1);
indElas=1;


p=nodes(:,2:3);
s=size(nodes,1);
%create adjacency matrix
adj=zeros(s,s);
adjSize=size(adj);
nodeLines=nodes(:,4:7);
for i=1:s
    for j=1:2
        if nodeLines(i,j)>0
            [row,~]=find(nodeLines(:,3:4)==nodeLines(i,j));
            adj(i,row)=1;
            %get linear index for fiber type
            if 0<nodeLines(i,j)<M_0
               if ismember(nodeLines(i,j),col)==1
                   collagen(indCol)= sub2ind(adjSize,i,row);
                   indCol=indCol+1;
               else
                   elastin(indElas)= sub2ind(adjSize,i,row);
                   indElas=indElas+1;
               end
            end
        end
    end
    for j=3:4
        if nodeLines(i,j)>0
            [row,~]=find(nodeLines(:,1:2)==nodeLines(i,j));
            adj(i,row)=1;
            %get linear index for fiber type
            if 0<nodeLines(i,j)<M_0
               if ismember(nodeLines(i,j),col)==1
                   collagen(indCol)= sub2ind(adjSize,i,row);
                   indCol=indCol+1;
               else
                   elastin(indElas)= sub2ind(adjSize,i,row);
                   indElas=indElas+1;
               end
            end
        end
    end
end

collagen=nonzeros(collagen);
elastin=nonzeros(elastin);

%% Graph
% G=graph(adj);
% [Ic, Jc]=ind2sub(size(adj), collagen);
% [Ie, Je]=ind2sub(size(adj), elastin);
% figure;
% h1=plot(G, 'XData',p(:,1),'YData',p(:,2));
% h1.EdgeAlpha=.8;
% highlight(h1,Ic,Jc,'EdgeColor',[0.6350, 0.0780, 0.1840], 'LineWidth', 2);
% highlight(h1,Ie,Je,'EdgeColor',[0, 0.4470, 0.7410], 'LineWidth', 1);
% pbaspect([1 1 1]);
% axis([0 1 0 1]);


end
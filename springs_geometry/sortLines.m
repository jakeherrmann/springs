function [rightLines, leftLines, topLines, bottomLines]=sortLines(xBeg, xEnd, yBeg, yEnd)
%% sortLines: find and sort which lines are attached to boundaries
%
%% INPUTS:
%
% xBeg == first x coordinate for each line
% xEnd == second x coordinate for each line
% yBeg == first y coordinate for each line
% yEnd == second y coordinate for each line
%
%% OUTPUTS
%
% for the following: ~Lines(1,2,3) = (value on boundary, line number, 0 or 1)
% where 1=beg and 0=end
% rightLines == line that have a node at x=1
% leftLines == line that have a node at x=0
% topLines == line that have a node at y=1
% bottomLines == line that have a node at y=0
%% LOGIC
%
% *this is this complicated because generateIsotropic Network needs all
% this to put springs on the boundaries*
%
% For each line, check if any have a node on the boundary
% Note its value, line number, and whether this is the beginning or end of
% the line
% Package up lines that belong to the same boundary and return
%
%% CODE
    bLineBeg=zeros(1); 
    tLineBeg=zeros(1);  
    rLineBeg=zeros(1);
    lLineBeg=zeros(1); 
    bLineEnd=zeros(1); 
    tLineEnd=zeros(1);  
    rLineEnd=zeros(1);
    lLineEnd=zeros(1);
%Sort lines
lBeg=1;
rBeg=1;
tBeg=1;
bBeg=1;
lEnd=1;
rEnd=1;
tEnd=1;
bEnd=1;

for s=1:size(xBeg,2)
   if xBeg(1,s) ==0
       lLineBeg(lBeg,2)=s;
       lLineBeg(lBeg,1)=yBeg(1,s);
       lLineBeg(lBeg,3)=1;
       lBeg=lBeg+1;
   end
   if xEnd(1,s)==0
       lLineEnd(lEnd,2)=s;
       lLineEnd(lEnd,1)=yEnd(1,s);
       lLineEnd(lEnd,3)=0;
       lEnd=lEnd+1;
   end
   if xBeg(1,s) ==1
       rLineBeg(rBeg,2)=s;
       rLineBeg(rBeg,1)=yBeg(1,s);
       rLineBeg(rBeg,3)=1;
       rBeg=rBeg+1;
   end
   if xEnd(1,s) ==1
       rLineEnd(rEnd,2)=s;
       rLineEnd(rEnd,1)=yEnd(1,s);
       rLineEnd(rEnd,3)=0;
       rEnd=rEnd+1;
   end
   if yBeg(1,s) ==0
       bLineBeg(bBeg,2)=s;
       bLineBeg(bBeg,1)=xBeg(1,s);
       bLineBeg(bBeg,3)=1;
       bBeg=bBeg+1;
   end
   if yEnd(1,s) ==0
       bLineEnd(bEnd,2)=s;
       bLineEnd(bEnd,1)=xEnd(1,s);
       bLineEnd(bEnd,3)=0;
       bEnd=bEnd+1;
   end
   if yBeg(1,s) ==1
       tLineBeg(tBeg,2)=s;
       tLineBeg(tBeg,1)=xBeg(1,s);
       tLineBeg(tBeg,3)=1;
       tBeg=tBeg+1;
   end
   if yEnd(1,s) ==1
       tLineEnd(tEnd,2)=s;
       tLineEnd(tEnd,1)=xEnd(1,s);
       tLineEnd(tEnd,3)=0;
       tEnd=tEnd+1;
   end
end
if size(lLineBeg,2)>1
    if size(lLineEnd,2)>1
        leftLines=[lLineBeg; lLineEnd];
    else
        leftLines=lLineBeg;
    end
else
    leftLines=lLineEnd;
end
if size(rLineBeg,2)>1
    if size(rLineEnd,2)>1
        rightLines=[rLineBeg; rLineEnd];
    else
        rightLines=rLineBeg;
    end
else
    rightLines=rLineEnd;
end
if size(tLineBeg,2)>1
    if size(tLineEnd,2)>1
        topLines=[tLineBeg; tLineEnd];
    else
        topLines=tLineBeg;
    end
else
    topLines=tLineEnd;
end
if size(bLineBeg,2)>1
    if size(bLineEnd,2)>1
        bottomLines=[bLineBeg; bLineEnd];
    else
        bottomLines=bLineBeg;
    end
else
    bottomLines=bLineEnd;
end
end
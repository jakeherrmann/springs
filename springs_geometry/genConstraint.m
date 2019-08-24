function [Aeq,beq]=genConstraint(n,edgeNodes)

Aeq = zeros(0, 2 * n);
beq = zeros(0, 1);

bot_idx=edgeNodes{3};
right_idx=edgeNodes{1};
top_idx=edgeNodes{2};
left_idx=edgeNodes{4};


for i = 1:length(right_idx)
    sz = size(Aeq, 1);
    Aeq(sz + 1, right_idx(i)) = 1;
    beq(sz+1) = 1;
end

for i = 1:length(top_idx)
    sz = size(Aeq, 1);
    Aeq(sz + 1, n + top_idx(i)) = 1;
    beq(sz+1) = 1;
end

for i = 1:length(bot_idx)
    sz = size(Aeq, 1);
    Aeq(sz + 1, n + bot_idx(i)) = 1;
    beq(sz+1) = 0;
end

for i = 1:length(left_idx)
    sz = size(Aeq, 1);
    Aeq(sz + 1, left_idx(i)) = 1;
    beq(sz+1) = 0;
end



beq = beq';

% Check
%check = all(Aeq * p(:) == beq);
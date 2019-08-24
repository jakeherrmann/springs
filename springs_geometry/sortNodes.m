function edgeNodes=sortNodes(p)

left_idx = find(p(:, 1) == 0);
right_idx = find(p(:, 1) == 1);
bot_idx = find(p(:, 2) == 0);
top_idx = find(p(:, 2) == 1);

edgeNodes{3}=bot_idx;
edgeNodes{1}=right_idx;
edgeNodes{2}=top_idx;
edgeNodes{4}=left_idx;
end

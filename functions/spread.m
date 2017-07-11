function Q = spread(grid, policy)

n_grid = length(grid);
Q = zeros(length(policy), n_grid);

for idx = 1:length(policy)
    left_idx = find(grid < policy(idx), 1, 'last');
    if (left_idx == n_grid)
        Q(idx,n_grid) = 1;
    elseif (isempty(left_idx))
        Q(idx,1) = 1;
    else
        to_left = (grid(left_idx+1) - policy(idx)) / (grid(left_idx+1) - grid(left_idx));
        Q(idx, left_idx) = to_left;
        Q(idx, left_idx+1) = 1-to_left;
    end
end
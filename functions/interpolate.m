function y = interpolate(x_grid, y_grid, x)
    left_idx = find(x_grid <= x, 1, 'last');
    if (left_idx == length(x_grid))
        y = myAD(x_grid(end), zeros(size(x.derivatives)));
    else
        slope = (y_grid(left_idx+1) - y_grid(left_idx))/(x_grid(left_idx+1) - x_grid(left_idx));
        y = y_grid(left_idx) + (x - x_grid(left_idx)) * slope;
    end
end
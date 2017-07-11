function plot_parameterwise(f, x, n, frac)

n_params = length(x);

for param_idx = 1:n_params
    
    x_grid = linspace(frac * x(param_idx), (1/frac) * x(param_idx), n);
    f_grid = zeros(n,1);
    for i = 1:n
        this_x = x;
        this_x(param_idx) = x_grid(i);
        f_grid(i) = f(this_x);
    end
    subplot(n_params, 1, param_idx);
    plot(x_grid, f_grid, '-o');
    yl = ylim();
    hold on;
    plot([x(param_idx) x(param_idx)], yl,'Color','red');
end
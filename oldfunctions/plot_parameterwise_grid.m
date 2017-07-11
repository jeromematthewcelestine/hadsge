function plot_parameterwise_grid(f, x, n)

n_params = length(x);

for param_idx = 1:n_params
    
    x_grid = linspace(0.01, 0.99, n);
    f_grid = zeros(n,1);
    for i = 1:n
        this_x = x;
        this_x(param_idx) = x_grid(i);
        f_grid(i) = f(this_x);
        f_grid(i)
    end
    subplot(n_params, 1, param_idx);
    plot(x_grid, f_grid, '-o');
%     yl = ylim();
%     hold on;
%     plot([x(param_idx) x(param_idx)], yl,'Color','red');
end
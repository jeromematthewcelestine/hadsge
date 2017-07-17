function kp_grid = solve_ss(params, k0, z0, Pz, kp_guess, options)

    nk = length(k0);
    nz = length(z0);

    policy_tolerance = options(1);
    kp = kp_guess;

    max_iter = 100;
    for iter = 1:max_iter
        kp_old = kp;

        max_abs_diff = 0;
        idx = 1;

        for z_idx = 1:nz
            Pz_row = Pz(:,z_idx);
            
            for k_idx = 1:nk

                % check lower bound
                foc_lb = foc(k0(1), params, z0, Pz_row, k0, k0(k_idx), kp_old);

                if (foc_lb > 0)
                    kp_star = k0(1);
                else
                    foc_ub = foc(k0(end), params, z0, Pz_row, k0, k0(k_idx), kp_old);

                    if (foc_ub < 0)
                        kp_star = k0(end);
                    else
                        kp_star = golden(k0(1), k0(end), params, options, z0, Pz_row, k0, k0(k_idx), kp_old);
                    end
                end

                kp(idx) = kp_star;

                abs_diff = abs(kp_star - kp_old(idx));
                if (abs_diff > max_abs_diff)
                    max_abs_diff = abs_diff;
                end

                idx = idx + 1;
            end
        end
        if (max_abs_diff < policy_tolerance)
            break;
        end
    end

    kp_grid = kp;
end

function x_star = golden(x_min, x_max, params, options, z0, Pz_row, k0, k, kp_grid)
    max_iter = 100;
    solver_tolerance = options(2);
    
    gr = 1.618;
    
    a = x_min;
    b = x_max;
    
    c = b - (b - a) / gr;
    d = a + (b - a) / gr;
    
    for i = 1:max_iter
        fc = abs(foc(c, params, z0, Pz_row, k0, k, kp_grid));
        fd = abs(foc(d, params, z0, Pz_row, k0, k, kp_grid));
        
        if (fc < fd)
            b = d;
        else
            a = c;
        end
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
        
        err = abs(c - d);
        
        if (err < solver_tolerance)
            break;
        end
    end
    x_star = (b + a) * 0.5;
end

function val = foc(kp, params, z0, Pz_row, k0, k, kp_grid)
    nz = length(z0);

    aalpha = params(1);
    bbeta = params(2);
    ddelta = params(3);
    pphi = params(4);
    
    one_minus_delta = 1 - ddelta;
    
    lhs = 1 + 2 * pphi * (kp/k - one_minus_delta);
    
    rhs = 0;
    for zp_idx = 1:nz
        kpp = interp(k0, kp_grid(:,zp_idx), kp);
        kpp_over_kp = kpp / kp;
        adj_part = pphi * (kpp_over_kp * kpp_over_kp - one_minus_delta * one_minus_delta);
        prodn_part = z0(zp_idx) * aalpha * kp^(aalpha-1.0);
        
        rhs = rhs + bbeta * Pz_row(zp_idx) * (prodn_part + one_minus_delta + adj_part);
    end
    val = lhs - rhs;
end

function yy = interp(x_grid, y_grid, xx)
    n = length(x_grid);
    
    if (xx <= x_grid(1))
        yy = y_grid(1);
        return;
    else
        for i = 2:n
            if (xx < x_grid(i))
                slope = (y_grid(i) - y_grid(i-1)) / (x_grid(i) - x_grid(i-1));
                yy = y_grid(i-1) + (xx - x_grid(i-1)) * slope;
                return;
            end
        end
        yy = y_grid(end);
        return;
    end
end

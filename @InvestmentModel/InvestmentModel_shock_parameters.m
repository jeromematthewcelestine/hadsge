classdef InvestmentModel < handle
properties
    % model parameters
    alpha % decreasing returns to scale
    beta % discount factor
    delta % depreciation rate
    phi % adjustment cost
    rho_z % persistence of idiosyncratic productivity
    sigma_z % volatility of idiosyncratic productivity
    rho_x_z % persistence of TFP shock
    sigma_x_z % volatility of TFP shock
    rho_x_b % persistence of discount factor shock
    sigma_x_b % volatility of discount factor shock
    rho_x_i % persistence of adjustment cost shock
    sigma_x_i % volatility of adjustment cost shock
    
    % misc
    idx 
    small_firm_index
    median_k
    
    n_states
    n_meas
    
    % z is idiosyncratic productivity, k is idiosyncratic capital
    nz
    nk
    z0
    Pz
    cum_Pz
    z_ergodic
    k0
    k_grid
    z_grid
    
    gensys_status
    
    % input matrices for Gensys
    G0
    G1
    Psi
    Pi
    Const
    
    % (solved) state space representation
    A
    B
    C
    D
    
    % steady state values
    Xss
end
methods
    
function o = InvestmentModel(params)
    o.alpha = params(1);
    o.beta = params(2);
    o.delta = params(3);
    o.phi = params(4);
    o.rho_z = params(5);
    o.sigma_z = params(6);
    o.rho_x_z = params(7);
    o.sigma_x_z = params(8);
    o.rho_x_b = params(9);
    o.sigma_x_b = params(10);
    o.rho_x_i = params(11);
    o.sigma_x_i = params(12);
    
    if (~o.check_bounds())
        o.gensys_status = [-3; -3];
        return;
    end 
    
    o.nk = 10;
    o.nz = 5;
    [o.Pz, o.z0] = rouwen(o.rho_z, 0.0, o.sigma_z, o.nz);
    o.Pz = o.Pz';
    o.cum_Pz = cumsum(o.Pz,2);
    o.z0 = exp(o.z0);
    
    T = 100;
    o.z_ergodic = (1/o.nz) * ones(o.nz,1);
    for t = 1:T
        o.z_ergodic = o.Pz' * o.z_ergodic;
    end
end

function solve(o)
    o.solve_steady_state();
    o.solve_derivatives();
    o.solve_dynamics();
end

function solve_steady_state(o)
    o.z_grid = reshape(repmat(o.z0', o.nk, 1), o.nk*o.nz, 1);
    options = [1e-6; 1e-3];

    k_max = 4;
    k_min = 1;

    max_n_tries = 5;
    for i = 1:max_n_tries;
        o.k0 = linspace(k_min, k_max, o.nk)';

        kp_guess = repmat(o.k0, 1, o.nz);
    
        o.k_grid = repmat(o.k0, o.nz, 1);

        kp_grid = solve_ss([o.alpha, o.beta, o.delta, o.phi], o.k0, o.z0, o.Pz', kp_guess, options);

        try_again = false;
        if (max(kp_grid(:)) >= o.k0(end))
            try_again = true;
            k_max = k_max * 2;
        elseif (max(kp_grid(:)) <= 0.5 * o.k0(end))
            try_again = true;
            k_max = k_max * 0.5;
        end
        if (min(kp_grid(:)) <= o.k0(1))
            try_again = true;
            k_min = k_min * 0.5;
        elseif (min(kp_grid(:)) >= 2.0 * o.k0(1))
            try_again = true;
            k_min = k_min * 2.0;
        end
        if (~try_again)
            break;
        end
    end
    
    kp_grid_tall = reshape(kp_grid, o.nk*o.nz, 1);

    Qk = spread(o.k0, kp_grid_tall);
    Q = repmat(Qk,1,o.nz) .* kron(o.Pz, ones(o.nk));

    dist = 1/(o.nk*o.nz) * ones(o.nk*o.nz,1);
    T = 100;
    for t = 1:T
        dist = Q' * dist;
    end

    median_k_idx = find(cumsum(sum(reshape(dist, o.nk, o.nz), 2)) > 0.25, 1);
    o.median_k = o.k0(median_k_idx);
    temp = reshape(1:o.nk*o.nz, o.nk, o.nz);
    o.small_firm_index = temp(1:median_k_idx,:);
    o.small_firm_index = o.small_firm_index(:);

    output_grid = o.z_grid .* o.k_grid .^ o.alpha;
    investment_grid = kp_grid_tall - (1-o.delta) .* o.k_grid;
    
    Yss = sum(output_grid(:) .* dist(:));
    Iss = sum(investment_grid(:) .* dist(:));
    I_small_ss = sum(investment_grid(o.small_firm_index) .* dist(o.small_firm_index));
    
    Css = Yss - Iss;
    pss = 1/Css;
    
    o.idx.x_z = 1;
    o.idx.x_b = 2;
    o.idx.x_i = 3;
    o.idx.kp = o.idx.x_i+1:o.idx.x_i+o.nz*o.nk;
    o.idx.Ekp = o.idx.kp(end)+1:o.idx.kp(end)+o.nk*o.nz;
    o.idx.dist = o.idx.Ekp(end)+1:o.idx.Ekp(end)+o.nk*o.nz;
    o.idx.p = o.idx.dist(end) + 1;
    o.idx.Ep = o.idx.p + 1;
    o.idx.Ex = o.idx.Ep + 1;
    o.idx.Ex_i = o.idx.Ex + 1;
    o.idx.y = o.idx.Ex_i + 1;
    o.idx.i = o.idx.y + 1;
    o.idx.c = o.idx.i + 1;
    o.idx.i_small = o.idx.c + 1;

    o.n_states = o.idx.i_small;

    o.Xss = zeros(o.n_states,1);

    o.Xss(o.idx.x_z) = 0;
    o.Xss(o.idx.x_b) = 0;
    o.Xss(o.idx.x_i) = 0;
    o.Xss(o.idx.kp) = kp_grid_tall;
    o.Xss(o.idx.Ekp) = kp_grid_tall;
    o.Xss(o.idx.dist) = dist;
    o.Xss(o.idx.p) = log(pss);
    o.Xss(o.idx.Ep) = log(pss);
    o.Xss(o.idx.Ex) = 0;
    o.Xss(o.idx.Ex_i) = 0;
    o.Xss(o.idx.y) = log(Yss);
    o.Xss(o.idx.i) = log(Iss);
    o.Xss(o.idx.c) = log(Css);
    o.Xss(o.idx.i_small) = log(I_small_ss);

    if (max(kp_grid(:)) >= o.k0(end))
        o.gensys_status = [-4; -4];
        return;
    end
end

function solve_derivatives(o)
    X_t = myAD([o.Xss; o.Xss]);
    
    X_t_r = reshape(X_t, o.n_states, 2);
    
    x_z_t   = exp(X_t_r(o.idx.x_z,1));
    x_b_t   = exp(X_t_r(o.idx.x_b,1));
    x_i_t   = exp(X_t_r(o.idx.x_i,1));
    kp_t    = X_t_r(o.idx.kp,1);
    Ekp_t   = X_t_r(o.idx.Ekp,1);
    dist_t  = X_t_r(o.idx.dist,1);
    p_t     = exp(X_t_r(o.idx.p,1));
    Ep_t    = exp(X_t_r(o.idx.Ep,1));
    Ex_t    = exp(X_t_r(o.idx.Ex,1));
    Ex_i_t    = exp(X_t_r(o.idx.Ex_i,1));
    Y_t     = exp(X_t_r(o.idx.y,1));
    I_t     = exp(X_t_r(o.idx.i,1));
    C_t     = exp(X_t_r(o.idx.c,1));
    I_small_t = exp(X_t_r(o.idx.i_small,1));
    
    x_z_tm1 = exp(X_t_r(o.idx.x_z,2));
    x_b_tm1 = exp(X_t_r(o.idx.x_b,2));
    x_i_tm1 = exp(X_t_r(o.idx.x_i,2));
    kp_tm1  = X_t_r(o.idx.kp,2);
    Ekp_tm1 = X_t_r(o.idx.Ekp,2);
    dist_tm1 = X_t_r(o.idx.dist,2);
    p_tm1   = exp(X_t_r(o.idx.p,2));
    Ep_tm1  = exp(X_t_r(o.idx.Ep,2));
    Ex_tm1  = exp(X_t_r(o.idx.Ex,2));
    Ex_i_tm1  = exp(X_t_r(o.idx.Ex_i,2));
    Y_tm1   = exp(X_t_r(o.idx.y,2));
    I_tm1   = exp(X_t_r(o.idx.i,2));
    C_tm1   = exp(X_t_r(o.idx.c,2));
    I_small_tm1 = exp(X_t_r(o.idx.i_small,2));
    
    o.Psi     = zeros(o.n_states,3);
    o.Psi(o.idx.x_z,1) = o.sigma_x_z;
    o.Psi(o.idx.x_b,2) = o.sigma_x_b;
    o.Psi(o.idx.x_i,3) = o.sigma_x_i;
    
    o.Pi      = zeros(o.n_states, o.nk * o.nz + 3);
    o.Pi(o.idx.Ekp, 1:o.nk * o.nz) = eye(o.nk * o.nz);
    o.Pi(o.idx.Ep, o.nk * o.nz +1) = 1;
    o.Pi(o.idx.Ex, o.nk * o.nz + 2) = 1;
    o.Pi(o.idx.Ex_i, o.nk * o.nz + 3) = 1;
    
    o.Const   = zeros(o.n_states,1);
    
    kp_grid_t_wide = reshape(kp_t, o.nk, o.nz);
    Ekp_grid_t_wide = reshape(Ekp_t, o.nk, o.nz);
    
    foc_wide = zeros(o.nk, o.nz) .* X_t(1);
    
    for k_idx = 1:o.nk
        k = o.k0(k_idx);
        for z_idx = 1:o.nz
            kp = kp_grid_t_wide(k_idx,z_idx);
            
            Erhs = zeros(o.nz,1) .* X_t(1);
            for zp_idx = 1:o.nz
                kpp = interpolate(o.k0, Ekp_grid_t_wide(:,zp_idx), kp);
                Erhs(zp_idx) = x_b_t * o.beta * ( Ex_t * o.z0(zp_idx) * o.alpha * kp.^(o.alpha-1) + 1 - o.delta + o.phi * Ex_i_t * ((kpp/kp).^2 - (1-o.delta)^2) );
            end
            rhs = Ep_t * sum(o.Pz(z_idx,:) .* Erhs');
            lhs = p_t * (1 + 2*o.phi*x_i_t*((kp./k) - (1-o.delta)));
            
            foc_wide(k_idx,z_idx) = lhs - rhs;
        end
    end
    foc = reshape(foc_wide, o.nk*o.nz, 1);
    
    Qk_new = zeros(o.nk*o.nz,o.nk) .* X_t(1);
    
    for i = 1:length(kp_t)
        left_idx = find(o.k0 < kp_t(i), 1, 'last');
        if (left_idx == o.nk)
            Qk_new(i, n_grid) = 1;
        elseif (isempty(left_idx))
            Qk_new(i, 1) = 1;
        else
            to_left = (o.k0(left_idx+1) - kp_t(i)) / (o.k0(left_idx+1) - o.k0(left_idx));
            Qk_new(i,left_idx) = to_left;
            Qk_new(i,left_idx+1) = 1-to_left;
        end
    end

    Q_new = repmat(Qk_new,1,o.nz) .* kron(o.Pz, ones(o.nk));
    
    dist_new = Q_new' * dist_tm1;
    
    output_grid_t = x_z_t * o.z_grid .* o.k_grid .^ o.alpha;
    investment_grid_t = kp_t - (1-o.delta) .* o.k_grid;
    
    Y_new = sum(dist_tm1(:) .* output_grid_t(:));
    I_new = sum(dist_tm1(:) .* investment_grid_t(:));
    I_small_new = sum(dist_tm1(o.small_firm_index) .* investment_grid_t(o.small_firm_index));
    C_new = Y_new - I_new;
    p_new = 1/C_new;

    f = zeros(o.n_states,1) * X_t(1);
    f(o.idx.x_z) = log(x_z_t) - o.rho_x_z*log(x_z_tm1);
    f(o.idx.x_b) = log(x_b_t) - o.rho_x_b*log(x_b_tm1);
    f(o.idx.x_i) = log(x_i_t) - o.rho_x_i*log(x_i_tm1);
    f(o.idx.kp) = foc;
	f(o.idx.Ekp) = Ekp_tm1 - kp_t;
	f(o.idx.dist)   = dist_t - dist_new;
	f(o.idx.p)      = log(p_t) - log(p_new);
	f(o.idx.Ep)     = log(Ep_tm1) - log(p_t);
    f(o.idx.Ex)     = log(Ex_tm1) - log(x_z_t);
    f(o.idx.Ex_i)     = log(Ex_i_tm1) - log(x_i_t);
	f(o.idx.y) = log(Y_t) - log(Y_new);
	f(o.idx.i) = log(I_t) - log(I_new);
	f(o.idx.c) = log(C_t) - log(C_new);
    f(o.idx.i_small) = log(I_small_t) - log(I_small_new);
    
    Gs = full(getderivs(f));
    o.G0 = Gs(:,1:o.n_states);
    o.G1 = -Gs(:,o.n_states+1:2*o.n_states);
    
    o.n_meas = 3;
    o.C = zeros(o.n_meas,o.n_states);
    o.C(1,o.idx.y) = 1;
    o.C(2,o.idx.i) = 1;
    o.C(3,o.idx.i_small) = 1;
end

function solve_dynamics(o)
    [o.A, ~, o.B, ~, ~, ~, ~, o.gensys_status] = gensys(o.G0, o.G1, o.Const, o.Psi, o.Pi);
end

function solve_dynamics_for_shock_parameters(o, shock_params)
    o.rho_x_z = shock_params(1);
    o.sigma_x_z = shock_params(2);
    o.rho_x_b = shock_params(3);
    o.sigma_x_b = shock_params(4);
    o.rho_x_i = shock_params(5);
    o.sigma_x_i = shock_params(6);
    
    o.Psi = zeros(o.n_states, 3);
    o.Psi(o.idx.x_z, 1) = o.sigma_x_z;
    o.Psi(o.idx.x_b, 2) = o.sigma_x_b;
    o.Psi(o.idx.x_i, 3) = o.sigma_x_i;
    o.G1(o.idx.x_z, o.idx.x_z) = o.rho_x_z;
    o.G1(o.idx.x_b, o.idx.x_b) = o.rho_x_b;
    o.G1(o.idx.x_i, o.idx.x_i) = o.rho_x_i;
    
    o.solve_dynamics();
end

function meas = simulate(o, innovations)
    T = size(innovations, 2);
    
    n_state = size(o.A, 1);
    state = zeros(n_state, T);
    
    o.n_meas = size(o.C, 1);
    meas = zeros(o.n_meas, T);
    
    meas(:, 1) = o.C * state(:, 1);
    for t = 2:T
        state(:, t) = o.A * state(:, t-1) + o.B * innovations(:, t);
        meas(:, t) = o.C * state(:, t); 
    end
    
    o.D = zeros(o.n_meas, o.n_meas);
end

% experimental; not yet in use
function [panel, meas] = simulate_panel(o, innovations, firm_innovations)
    T = size(innovations,2);
    
    n_state = size(o.A, 1);
    state = zeros(n_state,T);
    meas = zeros(o.n_meas,T);
    meas(:,1) = exp(o.C * o.Xss(:) + o.C * state(:,1));
    
    n_firms = size(firm_innovations,1);
    
    panel.k = zeros(n_firms,T);
    panel.z_idx = zeros(n_firms,T);
    
    panel.k(:,1) = 1;
    panel.z_idx(:,1) = 1;
    
    for t = 1:T
        if (t>1)
            state(:,t) = o.A * state(:,t-1) + o.B * innovations(:,t);
            meas(:,t) = exp(o.C * o.Xss(:) + o.C * state(:,t));
        end
        
        state_t = state(:,t);
        kp_t = reshape(o.Xss(o.idx.kp) + state_t(o.idx.kp), o.nk, o.nz);
        x_z_t = exp(o.Xss(o.idx.x_z));
        
        panel.I(t) = 0;
        panel.Is(t) = 0;
        panel.Y(t) = 0;
        
        for i = 1:n_firms
            z_idx = panel.z_idx(i,t);
            k = panel.k(i,t);
            kp = interpolate(o.k0, kp_t(:,z_idx), k);
            output = x_z_t * o.z0(z_idx) * k^o.alpha;
            investment = kp - (1-o.delta) * k;
            panel.investment(i,t) = investment;
            panel.I(t) = panel.I(t) + investment/n_firms;
            panel.Y(t) = panel.Y(t) + output/n_firms;
            if (k <= o.median_k)
                panel.Is(t) = panel.Is(t) + investment/n_firms;
            end
            panel.output(i,t) = output;
            if (t < T)
                panel.k(i,t+1) = kp;
                innov = firm_innovations(i,t+1);
                zp_idx = find(innov <= o.cum_Pz(z_idx,:), 1, 'first');
                panel.z_idx(i,t+1) = zp_idx;
            end
        end 
    end
end

function ll = evaluate_loglikelihood(o, data, innov_for_sim)    
    if (o.gensys_status(1)~=1)
        disp('gensys_status(1) ~= 1');
        ll = -Inf;
    else
        meas = data;
        ll = compute_loglikelihood(o.A, o.B, o.C, o.D, meas, innov_for_sim);
    end
end


function [logposterior, status] = evaluate_logposterior_for_parameters(o, params, data, innov_for_sim)
    phi         = params(1);
    rho_z       = params(2);
    sigma_z     = params(3);
    rho_x_z     = params(4);
    sigma_x_z   = params(5);
    rho_x_b     = params(6);
    sigma_x_b   = params(7);
    rho_x_i     = params(8);
    sigma_x_i   = params(9);

    % these priors are hard-coded for now
    logprior    = 0;
    logprior    = logprior + log(betapdf(rho_x_z,2,2));
    logprior    = logprior + log(inversegampdf(sigma_x_z,3,0.5));
    logprior    = logprior + log(betapdf(rho_x_b,2,2));
    logprior    = logprior + log(inversegampdf(sigma_x_b,3,0.5));
    logprior    = logprior + log(betapdf(rho_x_i,2,2));
    logprior    = logprior + log(inversegampdf(sigma_x_i,3,0.5));

    if (logprior > -Inf)
        [ll, status] = o.evaluate_loglikelihood_for_parameters(params, data, innov_for_sim);
        logposterior = logprior + ll;
    else
        status = -11;
        logposterior = -Inf;
    end
end

function distance = evaluate_smm_objective_function_for_shock_params(o, params, data_moments, innov_for_sim, invW)
        [sim_moments, status] = o.compute_simulated_moments_for_shock_params(params, innov_for_sim);
        if (status ~= 0)
            distance = Inf;
        else
            deviations = sim_moments - data_moments;
            distance = deviations(:)' * invW * deviations(:);
        end
    end

function [sim_moments,status] = compute_simulated_moments_for_shock_params(o, shock_params, innov_for_sim)
    o.solve_dynamics(shock_params);
    if (o.gensys_status(1) ~= 1)
        status = -1;
        sim_moments = [];
    else
        status = 0;
        sim_data = o.simulate(innov_for_sim);
        sim_moments = InvestmentModel.compute_moments(sim_data);
    end
end

function valid = check_bounds(o)
    valid = true;
    
    if (o.alpha <= 0 || o.alpha >= 1)
        valid = false;
        return;
    end
    if (o.beta <= 0 || o.beta >= 1)
        valid = false;
        return;
    end
    if (o.delta <= 0 || o.delta >= 1)
        valid = false;
        return;
    end
    if (o.phi < 0)
        valid = false;
        return;
    end
    if (o.rho_z <= 0 || o.rho_z >= 1)
        valid = false;
        return
    end
    if (o.sigma_z <= 0)
        valid = false;
        return
    end
    if (o.rho_x_z < 0 || o.rho_x_z >= 1)
        valid = false;
        return;
    end
    if (o.sigma_x_z < 0)
        valid = false;
        return;
    end
    if (o.rho_x_b < 0 || o.rho_x_b >= 1)
        valid = false;
        return;
    end
    if (o.sigma_x_b < 0)
        valid = false;
        return;
    end
    if (o.rho_x_i < 0 || o.rho_x_i >= 1)
        valid = false;
        return;
    end
    if (o.sigma_x_i < 0)
        valid = false;
        return;
    end
end

end

methods(Static)
    
function [ll, status] = evaluate_loglikelihood_for_parameters(params, data, innov_for_sim)
    new_model = InvestmentModel(params);
    status = sum(new_model.gensys_status);
    ll = new_model.evaluate_loglikelihood(data, innov_for_sim);
end

function [logposterior, status] = evaluate_logposterior_for_parameters(params, data, innov_for_sim)
    phi         = params(4);
    rho_z       = params(5);
    sigma_z     = params(6);
    rho_x_z     = params(7);
    sigma_x_z   = params(8);
    rho_x_b     = params(9);
    sigma_x_b   = params(10);
    rho_x_i     = params(11);
    sigma_x_i   = params(12);

    % these priors are hard-coded for now
    logprior    = log(gampdf(phi, 2, 1));
    logprior    = logprior + log(betapdf(rho_z,2,2));
    logprior    = logprior + log(inversegampdf(sigma_z,3,0.5));
    logprior    = logprior + log(betapdf(rho_x_z,2,2));
    logprior    = logprior + log(inversegampdf(sigma_x_z,3,0.5));
    logprior    = logprior + log(betapdf(rho_x_b,2,2));
    logprior    = logprior + log(inversegampdf(sigma_x_b,3,0.5));
    logprior    = logprior + log(betapdf(rho_x_i,2,2));
    logprior    = logprior + log(inversegampdf(sigma_x_i,3,0.5));

    if (logprior > -Inf)
        [ll, status] = InvestmentModel.evaluate_loglikelihood_for_parameters(params, data, innov_for_sim);
        logposterior = logprior + ll;
        fprintf('logprior %f, ll %f, logposterior %f\n', logprior, ll, logposterior);
    else
        status = -11;
        logposterior = -Inf;
    end
end


function [moments, W] = compute_moments(data)
    % compute moments for SMM
    data = data';

    mu = mean(data);

    n_rows = size(data,1);

    second_moments = (data - repmat(mu,n_rows,1)).^2;
    autocorrelations = (data(2:end,:) - repmat(mu,n_rows-1,1)) .* (data(1:end-1,:) - repmat(mu,n_rows-1,1));

    f = [second_moments(2:end,:) autocorrelations(1:end,:)];

    moments = mean(f,1);

    if (nargout == 2)
        n_rows = size(f,1);

        means = repmat(moments, n_rows, 1);
        wm = (f - means)';
        size(wm);
        W = (wm * wm') / (n_rows-1);
        maxq = 4;
        for iim = 1:maxq
            kapg = wm(:,iim+1:n_rows-1)*wm(:,1:n_rows-1-iim)';
            kapg = kapg + kapg';
            kapg = kapg/(n_rows-1);
            W = W + (1-iim/(maxq+1))*kapg;
        end
    end

    moments = moments';
end
end

end

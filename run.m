clear;
clc;

addpath('gensys');
addpath('functions');

alpha_star      = 0.3;
beta_star       = 0.9;
delta_star      = 0.05;
phi_star        = 0.1;
rho_z_star      = 0.9;
sigma_z_star    = 0.1;

rho_x_z_star    = 0.8;
sigma_x_z_star  = 0.2;
rho_x_b_star    = 0.8;
sigma_x_b_star  = 0.1;
rho_x_i_star    = 0.6;
sigma_x_i_star  = 0.1;

params_star = [alpha_star, beta_star, delta_star, phi_star, rho_z_star, sigma_z_star, rho_x_z_star, sigma_x_z_star, rho_x_b_star, sigma_x_b_star, rho_x_i_star, sigma_x_i_star];

% mex functions/solve_ss.c;

%% Setup and solve model

m = InvestmentModel(params_star);
m.solve();

%% Simulate model and plot

T = 200;
innov = randn(3,T);
data = m.simulate(innov);

plot(data')
l = legend('Y','I','Is');

%% Evaluate loglikelihood in neighborhood of true parameters

T_for_sim = 1000;
innov_for_sim = randn(3, T_for_sim);

posterior_obj_fn = @(theta) m.evaluate_logposterior_for_parameters(theta, data, innov_for_sim);
neg_obj_fn = @(theta) -posterior_obj_fn(theta);

figure;
plot_parameterwise(posterior_obj_fn, params_star(4:end), 5, 0.9)

%% Estimate using Metropolis-Hastings

chain_length = 10;

step_sizes = [0.01, 0.01, 0.01, 0.1, 0.01, 0.1, 0.01, 0.1, 0.01];
theta_init = [0.1, 0.9, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1];
param_names = {'rho x z', 'sigma x z', 'rho x b', 'sigma x b', 'rho x i', 'sigma x i'};

output_folder_name = 'mh_output';
if (~exist(output_folder_name,'dir'))
    mkdir(output_folder_name);
end

use_parallel_toolbox = true;

if (use_parallel_toolbox)
    
    n_inits = 2;
    theta_inits = zeros(n_inits,length(theta_init));
    theta_inits(1,:) = theta_init;
    theta_inits(2,:) = theta_init;

    spmd
        this_task = getCurrentTask();
        worker_id = this_task.ID;
        theta0 = theta_inits(mod(worker_id,n_inits)+1,:);

        file_id = fopen([output_folder_name,sprintf('/chain%04d.txt',worker_id)],'w');

        chain = metropolis_hastings(posterior_obj_fn, theta0, step_sizes, chain_length, file_id);

        fclose(file_id);
    end
else
    theta0 = theta_init;

    file_id = fopen([output_folder_name,sprintf('/chain%04d.txt',1)],'w');

    chain = metropolis_hastings(posterior_obj_fn, theta0, step_sizes, chain_length, file_id);

    fclose(file_id);
end


%% Load, combine and plot estimation results

results = MCMCResults('mh_output', 0.66);

n_params = 6;
for i = 1:n_params
    figure;
    hist(results.combined(:,i))
    xlabel(param_names{i})
end
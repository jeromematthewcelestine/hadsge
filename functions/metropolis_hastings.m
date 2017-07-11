function [chain, acceptance_rate, acceptance_rates] = metropolis_hastings(f, x0, step_sizes, chain_length, file_id)

theta_m1        = x0;
f_m1            = f(theta_m1);

n_accepts       = 0;

n_params = length(x0);

chain = zeros(chain_length,n_params);

param_counter = 1;

acceptances = zeros(1,n_params);
proposals = zeros(1,n_params);

for iter = 1:chain_length
    chain(iter,:)                 = theta_m1(:)';
    
    theta_star                  = theta_m1;
	theta_star(param_counter)   = theta_m1(param_counter) + step_sizes(param_counter) * randn;
    
    [f_star,status]             = f(theta_star);
    likelihood_ratio            = exp(f_star-f_m1);

    if (likelihood_ratio >= 1 || (f_m1 == -Inf && f_star > -Inf))
		accept = 1;
    elseif (likelihood_ratio < 0);
		accept = 0;
	else
		y = rand;
		if (y < likelihood_ratio)
			accept = 1;
		else
			accept = 0;
		end
    end
    
    proposals(param_counter) = proposals(param_counter) + 1;
    acceptances(param_counter) = acceptances(param_counter) + accept;
    
    if (accept == 1)
		theta_m1    = theta_star;
		f_m1        = f_star;
		n_accepts   = n_accepts + 1;
    end
    
    for i = 1:n_params
        fprintf(file_id, '%5.6f,', theta_m1(i));
    end
    fprintf(file_id, '%5.6f,%d\n',f_m1,accept);
    
    if (param_counter == n_params)
		param_counter = 1;
	else
		param_counter = param_counter + 1;
    end
end

acceptance_rates = acceptances ./ proposals;

acceptance_rate = n_accepts / chain_length;
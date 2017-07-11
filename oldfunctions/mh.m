function [chain, acceptance_rate] = mh(f, x0, step_sizes, chain_length)

theta_m1        = x0;
f_m1            = f(theta_m1);

n_accepts       = 0;

n_params = length(x0);

chain = zeros(chain_length,n_params);

param_counter = 1;

for iter = 1:chain_length
    chain(iter,:)                 = theta_m1(:)';
    
    theta_star                  = theta_m1;
% 	theta_star(param_counter)   = theta_m1(param_counter) + step_sizes(param_counter) * randn;
    theta_star                  = theta_m1 + step_sizes(:)' .* randn(1,n_params);
    f_star                      = f(theta_star);
    likelihood_ratio            = exp(f_star-f_m1);
	
    if (likelihood_ratio >= 1)
		accept = 1;
    elseif (likelihood_ratio < 0);
		accept = 0;
	else
		y = rand();
		if (y < likelihood_ratio)
			accept = 1;
		else
			accept = 0;
		end
    end

    if (accept == 1)
		theta_m1    = theta_star;
		f_m1        = f_star;
		n_accepts   = n_accepts + 1;
    end
    
    if (param_counter == n_params)
		param_counter = 1;
	else
		param_counter = param_counter + 1;
    end
end

acceptance_rate = n_accepts / chain_length;
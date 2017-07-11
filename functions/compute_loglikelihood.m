function ll = compute_loglikelihood(A, B, C, D, data, innov_for_sim)

n_state         = size(A,1);
data_length     = size(data, 2);

Sum_loglikelihood   = 0;
%__________________________________________________________________________
% Other matrices
Q               = B*B';
C               = C';
R				= D*D';
%__________________________________________________________________________
% Unconditional means and variance
P_1             = compute_unconditional_covariance_by_simulation(A, B, innov_for_sim);   %%%% Check
% P_1             = compute_unconditional_covariance_analytically(A, B)
% P_1             = compute_unconditional_covariance_by_simulation(A, B);   %%%% Check

xi_1            = zeros(n_state, 1);
%__________________________________________________________________________
% Initialize with x_10 and P_10 equal to unconditional mean and covariance
xi_10           = xi_1;     % unconditional mean of xi_1
P_10            = P_1;      % unconditional covariance
for t = 1:data_length
	%______________________________________________________________________
	mu                      = C'*xi_10;
	Sigma                   = C'*P_10*C + R;
	triuSigma               = triu(Sigma);
	Sigma                   = triuSigma + triuSigma' - diag(diag(Sigma));
    invSigma                = inv(Sigma);
	% DEBUG: Break if Sigma was not positive_definite
    [~,positive_definite]   = chol(Sigma);
	if (positive_definite~=0)
        fprintf('Error: Sigma not p.s.d at t = %d\n',t);
		disp(Sigma);
		ll = -Inf;
		return;
	end
	%______________________________________________________________________
    % Compute likelihood contribution
    Y_t                     = data(:,t);
	likelihood_t            = mvnpdf(Y_t', mu', Sigma);
	if (likelihood_t==0)
		Sum_loglikelihood = -Inf;
		break;
	end
	Sum_loglikelihood           = Sum_loglikelihood + log(likelihood_t);
	%______________________________________________________________________
    % Update mean forecasts
    xi_11                   = xi_10 + P_10*C*invSigma*(Y_t - C'*xi_10);
	xi_21                   = A * xi_11;                
	%______________________________________________________________________
    % Update variance of forecast of measurement variables
	P_11                    = P_10 - P_10*C*invSigma*C'*P_10;     
	P_21                    = A*P_11*A' + Q;                         
	%______________________________________________________________________
    % Update
	xi_10                   = xi_21;
	P_10                    = P_21;
end

ll = Sum_loglikelihood;
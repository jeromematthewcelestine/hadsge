function P = compute_unconditional_covariance_by_doubling(A, B)

Q = B*B';
S0  = Q;
d   = Inf;
tol = 1e-8;

while d>tol
    S1  = A*S0*A' + S0;
    A   = A*A;
	
    d   = max(max(abs(S1-S0)));
    err = max(max(abs(S0 - A*S0*A' - Q)));
    fprintf('d = %1.2e\t err = %1.2e\n',d,err);
    S0  = S1;
end

% OUTPUTS
P   = S1;
% err = max(max(abs(S - A*S*A' - Q)));
% fprintf('double error: %1.3e\n',err);
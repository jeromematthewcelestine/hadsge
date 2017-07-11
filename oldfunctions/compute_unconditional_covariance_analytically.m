function P = compute_unconditional_covariance_analytically(A, B)

tmp     = B*B';
tmp2    = (eye(size(A,1)^2) - kron(A,A)) \ tmp(:);
P       = reshape(tmp2,size(A,1),size(A,1));
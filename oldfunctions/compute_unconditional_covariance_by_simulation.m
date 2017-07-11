function P = compute_unconditional_covariance_by_simulation(varargin)

if (nargin >= 2)
    A = varargin{1};
    B = varargin{2};
end
if (nargin >= 3)
    innov = varargin{3};
    T = size(innov,2);
else
    T = 1000;
    innov = randn(size(B,2),T);
end

state = zeros(size(A,1),T);

state(:,1) = B * innov(:,1);
for t = 2:T
    state(:,t) = A * state(:,t-1) + B * innov(:,t);
end

P = cov(state');
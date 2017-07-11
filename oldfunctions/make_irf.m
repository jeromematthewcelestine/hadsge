function state = make_irf(A, B, shock_idx)

if (nargin==2)
	shock_idx=1;
end

T = 40;
state = zeros(size(A,1),T);
innov = zeros(size(B,2),T);
innov(shock_idx,2) = 1;

state(:,1) = B*innov(:,1);
for t = 2:T
	state(:,t) = A * state(:,t-1) + B * innov(:,t);
end
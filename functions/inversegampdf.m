function [ Y ] = inversegampdf( X,A,B )
%inversegampdf Inverse gamma probability density function.
%   Y = inversegampdf(X,A,B) returns the inverse gamma probability density
%   function with shape parameter A and scale parameter B, at the
%   values in X. The size of Y is the common size of the input arguments. A
%   scalar input functions is a constant matrix of the same size as the
%   other inputs.

Y = B.^A./gamma(A)*X.^(-A-1).*exp(-B./X);

end
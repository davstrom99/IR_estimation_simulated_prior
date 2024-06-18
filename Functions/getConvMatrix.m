function [X] = getConvMatrix(x,Nh)
% Takes an input signal and converts it into its equivalent matrix
% representation of the convolution operator.
%
% Input:
% x = Nx1, input signal
% Nh = Length of the impulse response used for the convolution operator
%
% Output:
% X = 

    N = length(x);    
    x_pad = [zeros(Nh-1,1); x; zeros(Nh-1,1)];
    X = zeros(N+Nh-1,Nh);

    for n = 1:N+Nh-1
        X(n,:) = flip(x_pad(n:n+Nh-1));
    end
end
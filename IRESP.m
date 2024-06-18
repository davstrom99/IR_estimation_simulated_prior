function [h_est] = IRESP(h_prior,X,y,eta,varargin)
%%% Impulse Response Estimation, Simulated Prior (IRESP)
%%% Computes the estimated room impulse response, using a simulated prior.
% [h_est] = SIOTBCD(h_prior,x,y,epsilon,eta,min_change_ratio, max_iter, max_iter_inner, C)
%
%%% Solves the problem 1/2 ||y - h*x||_2^2 + eta S(h,h_0), where * denotes
%%% the convolution operation, || . ||_2^2 is the 2-norm squared, and
%%% S(h,h_0) is the transport problem described in the paper (see end of
%%% function)
%
% Input: 
% 
% prior = Nh x 1 prior, a simulated RIR
% X = Convolution form of input signal, x (see getConvMatrix.m)
% y = (2*Nh+N-1) x 1 output-signal -- generated by conv(in-signal, IR)
% eta = regularization term
% epsilon = entropy term (default 1e-3)
% min_change_ratio = minimum ratio for the change of the estimated impulse
% response (default 0.005)
% max_iter = max number of iterations for the algorithm (default 120)
% max_iter_inner = max number of iterations for the inner loop (default 10)
% C = cost matrix, Nh x Nh (default based on euclidean distance)
%
% Output:
% 
% h_est = estimated impulse response, Nh x 1
%
%
% Based on the paper 'Room Impulse Response Estimation using Optimal
% Transport: Simulation-Informed Inference'
% https://arxiv.org/abs/2403.03762

    Nh = length(h_prior);

    % Initialize with default values
    epsilon = 1e-3;
    min_change_ratio = 0.005;
    max_iter = 120;
    max_iter_inner = 10;
    C = ((1:Nh)-(1:Nh)').^2;


    
    if nargin > 4
        if ~isempty(varargin{1}) && isnumeric(varargin{1}) && varargin{1} > 0
            epsilon = varargin{1};
        end
    end

    if nargin > 5
        if ~isempty(varargin{2}) && isnumeric(varargin{2}) && varargin{2} >= 0
            min_change_ratio = varargin{2};
        end
    end

    if nargin > 6
        if ~isempty(varargin{3}) && isnumeric(varargin{3}) && varargin{3} > 0
            max_iter = varargin{3};
        end
    end
    
    if nargin > 7
        if ~isempty(varargin{4}) && isnumeric(varargin{4}) && varargin{4} > 0
            max_iter_inner = varargin{4};
        end
    end

    if nargin > 8
        if ~isempty(varargin{5}) && all(all(varargin{5} >= 0))
            C = varargin{5};
        end
    end

    


    %Re-scaling the prior based on input and output (x, y)
    h_prior = scaleIRgeo(X,h_prior,y); 

    

    
    %step size
    gamma = 1/norm(X)^2;

    %initialize the estimated impulse response as zeros
    h_est = zeros(length(h_prior),1);

    hhat_old = h_est;
    mask_h_is = h_prior>0;

    dual0 = zeros(length(h_prior(mask_h_is)),1);
    dual1 = zeros(length(h_prior),1);

    tk = 1;
    rel_error = zeros(max_iter,1);
    for iter = 1:max_iter
        
        % Nestrov 'acceleration'
        tk_new = (1+sqrt(1+4*tk^2))/2;
        v = h_est + (tk-1)/tk_new*(h_est-hhat_old);
        tk = tk_new;

        % Gradient of LS term
        gradient = -X'*(y-X*v);

        [hhat_new,dual0,dual1] = prox_operator_square(h_prior,C,v-gamma*gradient,epsilon, eta, gamma,max_iter_inner,dual0,dual1);

        rel_error(iter) = norm(hhat_new-h_est)/norm(hhat_new);

        hhat_old = h_est;
        h_est = hhat_new;

        % Stopping criterion
        if iter > 1
            if rel_error(iter) < min_change_ratio
                return;
            end
        end

        if sum(isnan(h_est)) > 0 
            h_est = hhat_old;
            fprintf("Numerical instability \n")
            return
        end
        
    end
end

close all
clear all
clc

% files.mat contains two (room) impulse responses, both generated through
% the image source method. h_real is used to generate the signal, h_prior is
% used as a prior to estimate h_real. Also contains sampling frequency.
addpath Functions\
load("files.mat")
Nh = length(h_prior);
N = 350; %length of input signal

% input vector x
x = randn(N+Nh,1);

% recorded signal y, no noise
y_nn = conv(h_real,x);

% Adding some noise to the measured signal
signal_power = pow2db(mean(y_nn.^2));
y = awgn(y_nn, 10, signal_power);


% Taking the input signal and converting it into its equivalent matrix
% representation of the convolution operator, X.
[X] = getConvMatrix(x,Nh);
% Cutting the input and output signal
X = X(Nh+1:end-Nh+1,:);
x = x(Nh+1:end,:);
y = y(Nh+1:end-Nh+1,:);

figure(100)
hold off
plot(h_real)
hold on
plot(h_prior)
hold off
legend('Real RIR','Simulated prior RIR')


%% Estimation
eta = 1e-3;
h_est = IRESP(h_prior,X,y,eta);

%% Estimation with non-default values
% The values are chosen quite randomly and do not yield a good result; they
% are merely here to illustrate how to use the function.
% eta = 1e-7;
% epsilon = 1e-2;
% min_change_ratio = 0.0001;
% max_iter = 200;
% max_iter_inner = 25;
% C = (ones(1,Nh)'*(1:Nh));

% This takes much longer as the (maximum) number of iterations can be up to
% 4 times as many as for the default values
% h_est_example = IRESP(h_prior,X,y,eta, epsilon, min_change_ratio, max_iter, max_iter_inner, C);
% h_est_example = IRESP(h_prior,X,y,eta, [], [], max_iter, [], C);

%% Plotting
figure(101)
hold off
plot(h_real)
hold on
plot(h_est)
plot(h_prior)
legend('Real RIR', 'Estimated RIR','Simulated prior RIR')

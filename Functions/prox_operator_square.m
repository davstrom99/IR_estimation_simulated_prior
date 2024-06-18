function [h,dual0,dual1] = prox_operator_square(h_geo,C,u,epsilon, eta, gamma,max_iter,dual0,dual1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    nanvalues = false;
    theta = eta*gamma;

    mask_h_is = h_geo>0;

    rel_error_0 = zeros(max_iter,1);
    rel_error_1 = zeros(max_iter,1);
    
    ones_0 = ones(length(dual0),1);
    ones_1 = ones(length(dual1),1);

    dual0_new = zeros(size(dual0));
    dual1_new = zeros(size(dual1));

    iter = 1;
    while iter <= max_iter
        
        if nanvalues == false  %% if we don't need to use stable_log_sum_exp we don't
            % ---- dual0
            dual0_new = epsilon*theta*(log(h_geo(mask_h_is).^2) - log(exp(-1/epsilon*C(mask_h_is,:))*exp(1/(theta*epsilon)*dual1 ) ) );

            % ---- dual1
            dual1_new = 2*epsilon*theta*wrightOmegaq(1/2*log(u.^2) - ...
            1/2*log(16*(epsilon*theta)^2) - ...
            1/2* log(   (exp(-C(mask_h_is,:)/epsilon)).' * exp(dual0_new/(epsilon*theta)) ) + ...
            1/4/epsilon/theta) - 1/2;

            % max(...,0)
            dual1_new(dual1_new < 0) = 0;
        else %% if we have numerical instability we use stable_log_sum_exp
            % ---- dual0
            dual0_new = epsilon*theta*(log(h_geo(mask_h_is).^2) - stable_log_sum_exp( ...
                -C(mask_h_is,:)/epsilon + dual1.'/(epsilon*theta),ones_1));

            % ---- dual1
            dual1_new = 2*epsilon*theta*wrightOmegaq(1/2*log(u.^2) - ...
                1/2*log(16*(epsilon*theta)^2) - ...
                1/2*stable_log_sum_exp(-C(mask_h_is,:)'/epsilon + dual0_new'/(epsilon*theta),ones_0) ...
                + 1/4/epsilon/theta) - 1/2;

            % max(...,0)
            dual1_new(dual1_new < 0) = 0;
        end

        
        % Stopping criterion
        rel_error_0(iter) = norm(dual0-dual0_new)/norm(dual0_new);
        rel_error_1(iter) = norm(dual1-dual1_new)/norm(dual1_new);

        
        hasNaN = any(isnan( [dual0_new(:); dual1_new(:)] ));

        %---- We ignore this one iteration if there are NaN values
        if hasNaN == true
            nanvalues = true;
            dual0_new = dual0;
            dual1_new = dual1;
        end

        %saving the new dual variables
        dual1 = dual1_new;
        dual0 = dual0_new;

        h = u./(2*dual1+1);        

        % fprintf("Iteration %.d \n ",iter)
        if iter > 1
%             if rel_error_0(iter)<0.01 && rel_error_1(iter)<0.01
            if rel_error_0(iter)<0 && rel_error_1(iter)<0
                break;
            end
        end     
        if sum(isnan(dual1_new)) > 0
            h = h_geo*nan;
            fprintf('Numerical instability \n')
            return
        end
        
        if hasNaN == true
            iter = iter; %if we have nan values we do the iteration one more time and update the values with non-NaN values
        else
            iter = iter+1; %...otherwise continue as usual
        end
    end
    iter = iter-1;

    h = u./(2*dual1+1);
end
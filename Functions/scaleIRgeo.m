function [h_is] = scaleIRgeo(X,h_is,y)
% Scales the simulated impulse response such that it matches the energy in
% the recorded signal, y. 
    ir_scaling = (norm(y))/norm(X*h_is);
    h_is = h_is*ir_scaling;
end


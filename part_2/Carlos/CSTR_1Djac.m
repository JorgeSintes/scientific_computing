function J = CSTR_1Djac(t,x,p,F)
% Implementation of the CSTR 1D model
% x = Temperature in K

%we set the parameters
% p = [beta, k0, EaR, CAin, CBin, Tin, V];
%F = 0.02/60;
%if t>20
    F = F/60000;
    %end
beta = p(1); k0 = p(2); EaR = p(3); CAin = p(4); CBin = p(5); Tin = p(6); V = p(7);

J = (EaR*beta*k0*exp(-EaR/x)*(CAin + (Tin - x)/beta)*(CBin + (2*(Tin - x))/beta))/x^2 - 2*k0*exp(-EaR/x)*(CAin + (Tin - x)/beta) - k0*exp(-EaR/x)*(CBin + (2*(Tin - x))/beta) - F/V;


end
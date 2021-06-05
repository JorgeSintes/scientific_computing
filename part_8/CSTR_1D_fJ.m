function [Tdot,Tddot] = CSTR_1D_fJ(t, T, F, params)
% F to IS units
F = F/60000;

beta = params(1);
k0 = params(2);
EaR = params(3);
CAin = params(4);
CBin = params(5);
Tin = params(6);
V = params(7);

% Rate of reaction
r = k0 * exp(-EaR/T) * (CAin + 1/beta * (Tin - T)) * (CBin + 2/beta * (Tin - T));
% Rate of temperature
RT = beta * r;
% Final ODE
Tdot = F/V*(Tin - T) + RT;

% Derivative
Tddot = (EaR*beta*k0*exp(-EaR/T)*(CAin + (Tin - T)/beta)*(CBin + (2*(Tin - T))/beta))/T^2 ...
    - 2*k0*exp(-EaR/T)*(CAin + (Tin - T)/beta) - k0*exp(-EaR/T)*(CBin + (2*(Tin - T))/beta) - F/V;
end
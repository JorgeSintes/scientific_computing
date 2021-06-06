function [xdot, J] = CSTR_3D_fJ(t,x,F,params)
% F to IS units
F = F/60000;

beta = params(1);
k0 = params(2);
EaR = params(3);
CAin = params(4);
CBin = params(5);
Tin = params(6);
V = params(7);

CA = x(1);
CB = x(2);
T = x(3);
xdot = zeros(3,1);

% Rate of reaction
r = k0 * exp(-EaR/T) * CA * CB;
% Production rates and rate of temperature
RA = -r;
RB = -2*r;
RT = beta*r;
% Final ODE
xdot(1) = F/V * (CAin - CA) + RA;
xdot(2) = F/V * (CBin - CB) + RB;
xdot(3) = F/V * (Tin - T) + RT;

% Jacobian
J = zeros(3);
J(1,1) = -F/V - k0*exp(-EaR/T)*CB; J(1,2) = -k0*exp(-EaR/T)*CA; J(1,3) = -k0*EaR/T^2 *exp(-EaR/T)*CA*CB;
J(2,1) = -2*k0*exp(-EaR/T)*CB; J(2,2) = -F/V - 2*k0*exp(-EaR/T)*CA; J(2,3) = -2*k0*EaR/T^2 *exp(-EaR/T)*CA*CB;
J(3,1) = beta*k0*exp(-EaR/T)*CB; J(3,2) = beta*k0*exp(-EaR/T)*CA; J(3,3) = -F/V +beta*k0*EaR/T^2 *exp(-EaR/T)*CA*CB;
end
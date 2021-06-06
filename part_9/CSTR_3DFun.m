function [xdot,g] = CSTR_3DFun(t,x,args)
% F to IS units
[F, params] = args{:};
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

g = x;
end
function [xdot,J] = CSTR_3Dfunjac(t,x,p,F)
% Implementation of the CSTR 3D model
% x = [CA;CB;T]

%we set the parameters
% p = [beta, k0, EaR, CAin, CBin, Tin, V];
%F = 0.02/60;
%if t>20
    F = F/60000;
%end
xdot = zeros(3,1);

beta = p(1); k0 = p(2); EaR = p(3); CAin = p(4); CBin = p(5); Tin = p(6); V = p(7);
xdot(1) = F./V .* (CAin - x(1)) - k0.*exp(-EaR./x(3)).*x(1).*x(2); 
xdot(2) = F./V .* (CBin - x(2)) - 2.*k0.*exp(-EaR./x(3)).*x(1).*x(2); 
xdot(3) = F./V .* (Tin - x(3)) + beta*k0.*exp(-EaR./x(3)).*x(1).*x(2); 

J = zeros(3);
J(1,1) = -F/V - k0*exp(-EaR/x(3))*x(2); J(1,2) = -k0*exp(-EaR/x(3))*x(1); J(1,3) = -k0*EaR/x(3)^2 *exp(-EaR/x(3))*x(1)*x(2);
J(2,1) = -2*k0*exp(-EaR/x(3))*x(2); J(2,2) = -F/V - 2*k0*exp(-EaR/x(3))*x(1); J(2,3) = -2*k0*EaR/x(3)^2 *exp(-EaR/x(3))*x(1)*x(2);
J(3,1) = beta*k0*exp(-EaR/x(3))*x(2); J(3,2) = beta*k0*exp(-EaR/x(3))*x(1); J(3,3) = -F/V +beta*k0*EaR/x(3)^2 *exp(-EaR/x(3))*x(1)*x(2);
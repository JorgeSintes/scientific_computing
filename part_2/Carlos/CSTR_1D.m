function xdot = CSTR_1D(t,x,p,F)
% Implementation of the CSTR 1D model
% x = Temperature in K

%we set the parameters
% p = [beta, k0, EaR, CAin, CBin, Tin, V];
%F = 0.02/60;
%if t>20
    F = F/60000;
%end

beta = p(1); k0 = p(2); EaR = p(3); CAin = p(4); CBin = p(5); Tin = p(6); V = p(7);
xdot = F./V .* (Tin - x) + beta.*k0.*exp(-EaR./x).*(CAin +(1./beta).*(Tin-x)) .*(CBin +(2./beta).*(Tin-x));
end
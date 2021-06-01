function g = CSTR_3D_Diffusion(t,x,F,params,sigma)
% F to IS units
F = F/60000;
V = params(7);

g = zeros(3,1);

g(3) = F/V * sigma;
end
function g = CSTR_1D_Diffusion(t,x,F,params,sigma)
% F to IS units
F = F/60000;
V = params(7);

g = F/V * sigma;
end
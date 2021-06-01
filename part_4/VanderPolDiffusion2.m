function g = VanderPolDiffusion2(t,x,mu,sigma)
g = [0; sigma*(1.0+x(1)^2)];
end
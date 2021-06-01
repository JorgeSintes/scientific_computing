function [Tw,W,dW] = StdWienerProcess(tspan,h,nW,Ns,seed)

if nargin == 5
    rng(seed);
end

t0 = tspan(1);
tf = tspan(2);

Tw = t0:h:tf;
N = size(Tw,2);

dW = sqrt(h)*randn(nW,N,Ns);
W  = [zeros(nW,1,Ns) cumsum(dW,2)];
W = W(:,1:end-1,:);
end

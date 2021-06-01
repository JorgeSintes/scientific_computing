function [T,X] = ImplicitExplicit(ffun,gfun,tspan,h,x0,nW,args)

newtontol = 1e-8;
maxit = 100;

t0 = tspan(1);
tf = tspan(end);
T = t0:h:tf;
N  = size(T,2);

X  = zeros(size(x0,1),N);
X(:,1) = x0;

[~,~,dW] = StdWienerProcess(tspan,h,nW,1);
k = 1;

f = feval(ffun,T(k),X(:,k),args{:});
for k=1:N-1
    g = feval(gfun,T(k),X(:,k),args{:});

    psi = X(:,k) + g*dW(:,k);

    xinit = psi + h*f;
    [X(:,k+1),f,~] = NewtonsMethodSDE(ffun,T(:,k+1),h,psi,xinit,newtontol,maxit,args);
end
end
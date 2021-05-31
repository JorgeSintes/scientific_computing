function [T,X] = EulerImplicit_fixed(funJac, tspan, h, x0, args)

t0 = tspan(1);
tf = tspan(end);
T = t0:h:tf;
N = size(T,2);
X = zeros(size(x0,1), N);
X(:,1) = x0;

tol = 1.0e-8;
maxit = 100;

for k=1:N-1
    f = feval(funJac,T(k),X(:,k),args{:});
    xinit = X(:,k) + f*h;
    [X(:,k+1),~] = NewtonsMethod(funJac, T(:,k), X(:,k),...
        h, xinit, tol, maxit, args);
end
end

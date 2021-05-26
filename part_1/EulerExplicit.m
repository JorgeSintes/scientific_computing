function [T,X,X_real,e_l,e_g] = EulerExplicit(fun, anal_sol, tspan, h, x0, args)

t0 = tspan(1);
tf = tspan(end);
T = t0:h:tf;
N = size(T,2);
X = zeros(size(x0,1), N);
X_real = zeros(size(x0,1), N);
e_g = zeros(size(X));
e_l = zeros(size(X));
X(:,1) = x0;
X_real(:,1) = x0;

for k = 1:N-1
    f = feval(fun, T(k), X(:,k), args{:});
    X(:,k+1) = X(:,k) + h * f;
    X_real(:,k+1) = feval(anal_sol, T(k+1), X(:,1), T(1), args{:});
    
    e_l(:,k+1) = abs(X(:,k+1) - feval(anal_sol, T(k+1), X(:,k), T(k), args{:}));
    e_g(:,k+1) = abs(X(:,k+1) - X_real(:,k+1));
end
end
function [T,X] = EulerExplicit_fixed(fun, tspan, h, x0, args)

t0 = tspan(1);
tf = tspan(end);
T = t0:h:tf;
N = size(T,2);
X = zeros(size(x0,1), N);
X(:,1) = x0;

for k = 1:N-1
    f = feval(fun, T(k), X(:,k), args{:});
    X(:,k+1) = X(:,k) + h * f;
end
end

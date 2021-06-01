function [T,X] = EulerMaruyama(ffun,gfun,tspan,h,x0,nW,args)
t0 = tspan(1);
tf = tspan(end);
T = t0:h:tf;
N  = size(T,2);

X  = zeros(size(x0,1),N);
X(:,1) = x0;

[~,~,dW] = StdWienerProcess(tspan,h,nW,1);

for k=1:N-1
    f = feval(ffun,T(k),X(:,k),args{:});
    g = feval(gfun,T(k),X(:,k),args{:});

    psi = X(:,k) + g*dW(:,k);

    X(:,k+1) = psi + h*f;
end

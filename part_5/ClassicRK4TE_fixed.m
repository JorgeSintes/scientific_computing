function [T,X,X_real,e_l,e_g] = ClassicRK4TE_fixed(fun,anal_sol,tspan,h,x0,args)
% time interval
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
    t = T(k);
    x = X(:,k);

    % Stage 1
    T1 = t;
    X1 = x;
    F1 = feval(fun,T1,X1,args{:});
    % Stage 2
    T2 = t + h/2;
    X2 = x + h/2*F1;
    F2 = feval(fun,T2,X2,args{:});
    % Stage 3
    T3 = T2;
    X3 = x + h/2*F2;
    F3 = feval(fun,T3,X3,args{:});
    % Stage 4
    T4 = t + h;
    X4 = x + h*F3;
    F4 = feval(fun,T4,X4,args{:});

    % final solution
    X(:,k+1) = x + h*(1/6*F1+1/3*F2+1/3*F3+1/6*F4);

    X_real(:,k+1) = feval(anal_sol, T(k+1), x0, t0, args{:});
    
    e_l(:,k+1) = abs(X(:,k+1) - feval(anal_sol, T(k+1), X(:,k), T(k), args{:}));
    e_g(:,k+1) = abs(X(:,k+1) - X_real(:,k+1));
   
end
end
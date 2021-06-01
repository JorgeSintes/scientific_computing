function xout = ClassicRK4_step(fun,t,h,x,args)
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
xout = x + h*(1/6*F1+1/3*F2+1/3*F3+1/6*F4);
end

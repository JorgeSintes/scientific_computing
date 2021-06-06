function [x,k] = NewtonsMethod(FunJac, tk, xk, h, xinit, tol, maxit, args)
k = 0;
t = tk + h;
x = xinit;
[f,J] = feval(FunJac,t,x,args{:});
k = k + 1;
R = x - h*f - xk;
I = eye(length(xk));

while( (k<maxit) && (norm(R,'inf')>tol) )
    dRdx = I - J*h;
    dx = dRdx\R;
    x = x - dx;
    [f,J] = feval(FunJac,t,x,args{:});
    k = k+1;
    R = x - h*f - xk;
end

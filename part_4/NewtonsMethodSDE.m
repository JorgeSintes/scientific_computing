function [x,f,J] = NewtonsMethodSDE(ffun,t,h,psi,xinit,tol,maxit,args)
I = eye(length(xinit));
x = xinit;

[f,J] = feval(ffun,t,x,args{:});
R = x - h*f - psi;
it = 1;

while ( (norm(R,'inf')>tol) &&(it<= maxit) )
    dRdx = I - J*h;
    mdx = dRdx\R;

    x   = x - mdx;
    [f,J] = feval(ffun,t,x,args{:});
    R = x - h*f - psi;
    it = it + 1;
end
end
function [T,X,r_out,h_out,info] = EulerImplicit_adaptive(funJac,tspan,h0,x0,abstol,reltol,args)

epstol = 0.8;
facmin = 0.1;
facmax = 1.5;

newtontol = 1.0e-8;
maxit = 100;

t0 = tspan(1);
tf = tspan(end);
t = t0;
h = h0;
x = x0;

T = t0;
X = x0;
r_out = [];
h_out = [];
info = zeros(1,4);

nfun = 0;
nstep = 0;
naccept = 0;

while t < tf
    if (t+h > tf)
        h = tf-t;
    end
    [f,J] = feval(funJac,t,x,args{:});
    nfun = nfun + 1;
    
    AcceptStep = false;
    while ~AcceptStep
        % Single step size
        xinit1 = x + h*f;
        [x1,nfun_local] = NewtonsMethod(funJac, t, x, h, xinit1, newtontol, maxit, args);
        nfun = nfun + nfun_local;
        
        hm = 0.5*h;
        tm = t + hm;
        xinitm = x + hm*f;
        [xm,nfun_local] = NewtonsMethod(funJac, t, x, hm, xinitm, newtontol, maxit, args);
        nfun = nfun + nfun_local;
        
        [fm,Jm] = feval(funJac,tm,xm,args{:});
        nfun = nfun + 1;
        xinit1hat = xm + hm*f;
        [x1hat,nfun_local] = NewtonsMethod(funJac, tm, xm, hm, xinit1hat, newtontol, maxit, args);
        nfun = nfun + nfun_local;
        
        nstep = nstep + 1;
        
        e = abs(x1hat-x1);
        r = max(e./max(abstol, abs(x1hat) .* reltol));
        AcceptStep = (r <= 1.0);
        
        if AcceptStep
            t = t+h;
            x = x1hat;
           
            T = [T,t];
            X = [X,x];
            r_out = [r_out, r];
            h_out = [h_out, h];
            
            naccept = naccept + 1;
        end
        
        h = max(facmin, min(sqrt(epstol/r), facmax)) * h;
    end
end

info(1) = nfun;
info(2) = nstep;
info(3) = naccept;
info(4) = nstep - naccept;

end
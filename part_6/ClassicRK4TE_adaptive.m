function [T,X,X_real,e_l,e_g,r_out,h_out,info] = ClassicRK4TE_adaptive(fun,anal_sol,tspan,h0,x0,abstol,reltol,args)

epstol = 0.8;
kpow = 0.2;
facmin = 0.1;
facmax = 5.0;

t0 = tspan(1);
tf = tspan(end);
t = t0;
h = h0;
x = x0;

T = t0;
X = x0;
X_real = x0;
e_l = 0;
e_g = 0;
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
    
    nfun = nfun + 1;
    AcceptStep = false;
    while ~AcceptStep
        x1 = ClassicRK4_step(fun,t,h,x,args);
        nfun = nfun + 4;
        
        hm = 0.5*h;
        tm = t + hm;
        xm = ClassicRK4_step(fun,t,hm,x,args);
        nfun = nfun + 4;
        x1hat = ClassicRK4_step(fun,tm,hm,xm,args);
        nfun = nfun + 4;
        nstep = nstep + 1;
        
        e = abs(x1hat-x1);
        r = max(e./max(abstol, abs(x1hat) .* reltol));
        AcceptStep = (r <= 1.0);
        
        if AcceptStep
            t = t+h;
            x = x1hat;
            
            x_real = feval(anal_sol, t, x0, t0, args{:});
                        
            e_l = [e_l, abs(x - feval(anal_sol, t, X(:,end), T(end), args{:}))];
            e_g = [e_g, abs(x - x_real)];
            
            T = [T,t];
            X = [X,x];
            X_real = [X_real, x];
            r_out = [r_out, r];
            h_out = [h_out, h];
            naccept = naccept + 1;
        end
        
        h = max(facmin, min((epstol/r)^kpow, facmax)) * h;
    end
end

info(1) = nfun;
info(2) = nstep;
info(3) = naccept;
info(4) = nstep - naccept;

end
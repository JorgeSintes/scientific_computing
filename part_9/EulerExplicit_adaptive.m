function [T,X,r_out,h_out,info] = EulerExplicit_adaptive(fun,tspan,h0,x0,abstol,reltol,args)

epstol = 0.8;
facmin = 0.1;
facmax = 5.0;

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
    f = feval(fun,t,x,args{:});
    nfun = nfun + 1;
    AcceptStep = false;
    while ~AcceptStep
        x1 = x + h*f;
        
        hm = 0.5*h;
        tm = t + hm;
        xm = x + hm*f;
        fm = feval(fun,tm,xm,args{:});
        nfun = nfun + 1;
        x1hat = xm + hm*fm;
        
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

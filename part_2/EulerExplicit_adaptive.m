function [T,X] = EulerExplicit_adaptive(fun,tspan,h0,x0,abstol,reltol,args)

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

while t < tf
    if (t+h > tf)
        h = tf-t;
    end
    f = feval(fun,t,x,args{:});
    
    AcceptStep = false;
    while ~AcceptStep
        x1 = x + h*f;
        
        hm = 0.5*h;
        tm = t + hm;
        xm = x + hm*f;
        fm = feval(fun,tm,xm,args{:});
        x1hat = xm + hm*fm;
        
        e = abs(x1hat-x1);
        r = max(e./max(abstol, abs(x1hat) .* reltol));
        AcceptStep = (r <= 1.0);
        
        if AcceptStep
            t = t+h;
            x = x1hat;
            
            T = [T,t];
            X = [X,x];
        end
        
        h = max(facmin, min(sqrt(epstol/r), facmax)) * h;
    end
end

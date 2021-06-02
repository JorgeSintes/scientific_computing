function [T_out,X_out,e_out,r_out,h_out,info] = ... 
    DOPRI54_adaptive(fun,tspan,h0,x0,abstol,reltol,butcher,args)

% Controller parameters
epstol = 0.8;
kpow = 1/6;
facmin = 0.1;
facmax = 5.0;

t0 = tspan(1);
tf = tspan(end);
t = t0;
h = h0;
x = x0;

nx = size(x0,1);

T_out = t0;
X_out = x0;
X_real = x0;
e_out = zeros(nx,1);
r_out = [];
h_out = [];
info = zeros(4,1);

nfun = 0;
nstep = 0;
naccept = 0;

% Extract Butcher Tableau
s = butcher.stages;
AT = butcher.AT;
b = butcher.b;
c = butcher.c;
d = butcher.d;

% Allocate memory for the stages
T = zeros(1,s);
X = zeros(nx,s);
F = zeros(nx,s);

while t < tf
    
    % Size of last step
    if t+h > tf
        h = tf - t;
    end
    % First stage
    T(1) = t;
    X(:,1) = x;
    F(:,1) = feval(fun,T(1),X(:,1),args{:});
    nfun = nfun + 1;
    
    AcceptStep = false;
   
    while ~AcceptStep
        % Precalculated parameters
        hAT = h*AT;
        hb = h*b;
        hc = h*c;
        hd = h*d;

        % Following stages
        for i = 2:s
            T(i) = t + hc(i);
            X(:,i) = x + F(:,1:i-1)*hAT(1:i-1,i);
            F(:,i) = feval(fun,T(i),X(:,i),args{:});
            nfun = nfun + 1;
        end
        
        % Error estimation
        e = F*hd;
        xhat = x + F*hb;
        r = max(abs(e)./max(abstol,abs(xhat).*reltol));
         % Check conditon
        AcceptStep = r <=1;
        
        if AcceptStep
            t = t+h;
            x = xhat;
            
            T_out = [T_out,t];
            X_out = [X_out,x];
            e_out = [e_out,e];
            r_out = [r_out, r];
            h_out = [h_out, h];
            naccept = naccept + 1;
        end
        
        % step size controller (Asymptotic or second order PI)
        h = max(facmin,min((epstol/r)^kpow,facmax))*h;
        nstep = nstep+1;    
    end
end

info(1) = nfun;
info(2) = nstep;
info(3) = naccept;
info(4) = nstep - naccept;

end
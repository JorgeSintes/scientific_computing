function [T,X,r_eul,h_eul,info] = ImplicitEulerAdaptiveStepSize(funJac,tspan,x0,h0,abstol,reltol,varargin)
% Error controller parameters
epstol = 0.8; % target
facmin = 0.1; % maximum decrease factor
facmax = 5.0; % maximum increase factor
tol = 1.0e-8;
maxit = 100;
% Integration interval
t0 = tspan(1);
tf = tspan(2);
% Initial conditions
t = t0;
h = h0;
x = x0;
% Output
T = t;
X = x';
%% Main algorithm
r_eul = [0];
h_eul = [h0];
nfun = 0;
accepted = 0;
while t < tf
    if (t+h>tf)
        h = tf-t;
    end
    f = feval(funJac,t,x,varargin{:});
    nfun = nfun + 1;
    
    AcceptStep = false;
    while ~AcceptStep
        % Take step of size h
        xinit1 = x + h*f;
        [x1,nnfun] = NewtonsMethodODE(funJac,...
        t, x, h, xinit1, tol, maxit, varargin{:});  
        nfun = nfun + nnfun;
        
        % Take step of size h/2
        hm = 0.5*h;
        tm = t + hm;
        xinitm = x + hm*f;
        [xm,nnfun] = NewtonsMethodODE(funJac,...
        t, x, hm, xinitm, tol, maxit, varargin{:});  
        nfun = nfun + nnfun;
        
        fm = feval(funJac,tm,xm,varargin{:});
        nfun = nfun + 1;
        
        xinit1hat = xm + hm*fm;
        [x1hat,nnfun] = NewtonsMethodODE(funJac,...
        tm, xm, hm, xinit1hat, tol, maxit, varargin{:});  
        nfun = nfun + nnfun;
        
        % Error estimation
        e = x1hat-x1;
        r = max(abs(e)./max(abstol,abs(x1hat).*reltol));
        AcceptStep = (r<=1.0);
        if AcceptStep
            t = t+h;
            x = x1hat;
            
            T = [T;t];
            X = [X;x'];
            h_eul = [h_eul;max(facmin,min(sqrt(epstol/r),facmax))*h];
            r_eul = [r_eul;r];
            accepted = accepted + 1;
        end
        %Asymptotic step size controller
        h = max(facmin,min(sqrt(epstol/r),facmax))*h;

    end

    info.nfun = nfun;
    info.accepted = accepted;
    
end


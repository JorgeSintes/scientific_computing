%% CSTR 3D implicit for different tolerances
clc
clear
close all

% Parameter setting
beta = 560/(1.0*4.186);
k0 = exp(24.6);
EaR = 8500;
CAin = 1.6/2;
CBin = 2.4/2;
Tin = 273.65;
V = 0.105;

% Piecewise flow definition
tspans = [[0,3];[3,5];[5,7];[7,9];[9,12];[12,16];[16,18];...
    [18,20];[20,22];[22,24];[24,28];[28,32];[32,35]];
Fs = [700,600,500,400,300,200,300,400,500,600,700,200,700];

% Flow for the plot
tspans_init = [];
F_init = [];
for i=1:size(Fs,2)
    tspans_init = [tspans_init,tspans(i,:)];
    F_init = [F_init,Fs(i), Fs(i)];   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [1e-2 1e-3 1e-4 1e-5 1e-6];

count = 0;
infos = zeros(4,size(tols,2));

figure('Renderer', 'painters', 'Position', [10 10 1000 800])

for tol=tols
    count = count+1;
    
    x0 = [[0; 0; Tin]];
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,r_local,h_local, info_local] = EulerExplicit_adaptive(@CSTR_3D,tspan,h0,x0,abstol,reltol,args);
        
        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];
        r_3D = [r_3D, r_local(:,1:end)];
        h_3D = [h_3D, h_local(:,1:end)];
        info_3D = info_3D + info_local';
        
        x0 = X_local(:,end);
        h0 = h_local(end);
    end
    
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    infos(:,count) = info_3D;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on
xlabel('tol');
ylabel('Function Evaluations');
set(gca,'Xdir','reverse');

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on
xlabel('tol');
ylabel('Calculated steps');
set(gca,'Xdir','reverse');

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on
xlabel('tol');
ylabel('Accepted steps');
set(gca,'Xdir','reverse');

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on
xlabel('tol');
ylabel('Rejected steps');
set(gca,'Xdir','reverse');

%% Implicit

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count+1;
    
    x0 = [[0; 0; Tin]];
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,r_local,h_local, info_local] = EulerImplicit_adaptive(@CSTR_3D_fJ,tspan,h0,x0,abstol,reltol,args);
        
        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];
        r_3D = [r_3D, r_local(:,1:end)];
        h_3D = [h_3D, h_local(:,1:end)];
        info_3D = info_3D + info_local';
        
        x0 = X_local(:,end);
        h0 = h_local(end);
    end
    
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    infos(:,count) = info_3D;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% Classical Runge-Kutta


count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count+1;
    
    x0 = [[0; 0; Tin]];
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,r_local,h_local, info_local] = ClassicRK4_adaptive(@CSTR_3D_fJ,tspan,h0,x0,abstol,reltol,args);
        
        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];
        r_3D = [r_3D, r_local(:,1:end)];
        h_3D = [h_3D, h_local(:,1:end)];
        info_3D = info_3D + info_local';
        
        x0 = X_local(:,end);
        h0 = h_local(end);
    end
    
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    infos(:,count) = info_3D;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% DOPRI54

butcher = ERKSolverErrorEstimationParameters('DOPRI54');
count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count+1;
    
    x0 = [[0; 0; Tin]];
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,e,r_local,h_local, info_local] = DOPRI54_adaptive(@CSTR_3D,tspan,h0,x0,abstol,reltol,butcher,args);
        
        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];
        r_3D = [r_3D, r_local(:,1:end)];
        h_3D = [h_3D, h_local(:,1:end)];
        info_3D = info_3D + info_local;
        
        x0 = X_local(:,end);
        h0 = h_local(end);
    end
    
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    infos(:,count) = info_3D;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% ESDIRK23

method = 'ESDIRK23';

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count+1;
    
    x0 = [[0; 0; Tin]];
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        t0 = tspan(1);
        tf = tspan(end);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,e,struct_info_local,stats_local] = ESDIRK(@CSTR_3DFun,@CSTR_3DJac,tspan(1),tspan(2),x0,h0,abstol,reltol,method,args);
        X_local = X_local';
        T_local = T_local';
        
        info_local(1) = struct_info_local.nFun;
        info_local(2) = struct_info_local.nStep;
        info_local(3) = struct_info_local.nAccept;
        info_local(4) = struct_info_local.nFail;
        
        info_3D = info_3D + info_local;
        
        x0 = X_local(:,end);
        h0 = h_local(end);
    end
    
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    infos(:,count) = info_3D;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% ode45

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count+1;
    
    x0 = [[0; 0; Tin]];
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        t0 = tspan(1);
        tf = tspan(end);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        
        opts = odeset('RelTol',tol,'AbsTol',tol);
        sol_3D_45 = ode45(@(t,x) CSTR_3D(t,x,F,params),tspan, x0, opts);
        
        info_local(1) = info_local(1) + sol_3D_45.stats.nfevals;
        info_local(3) = info_local(3) + sol_3D_45.stats.nsteps;
        info_local(4) = info_local(4) + sol_3D_45.stats.nfailed;
        info_local(2) = info_local(2) + info_local(3) + info_local(4);
        
        info_3D = info_3D + info_local;
        
        x0_45 = sol_3D_45.y(:,end);
    end
    
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    infos(:,count) = info_3D;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% ode15s

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count+1;
    
    x0 = [[0; 0; Tin]];
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        t0 = tspan(1);
        tf = tspan(end);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        
        opts = odeset('RelTol',tol,'AbsTol',tol);
        sol_3D_45 = ode15s(@(t,x) CSTR_3D(t,x,F,params),tspan, x0, opts);
        
        info_local(1) = info_local(1) + sol_3D_45.stats.nfevals;
        info_local(3) = info_local(3) + sol_3D_45.stats.nsteps;
        info_local(4) = info_local(4) + sol_3D_45.stats.nfailed;
        info_local(2) = info_local(2) + info_local(3) + info_local(4);
        
        info_3D = info_3D + info_local;
        
        x0_45 = sol_3D_45.y(:,end);
    end
    
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    infos(:,count) = info_3D;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on
legend('Ex. Euler', 'Im. Euler','RK4','DOPRI54','ESDIRK23','ode45','ode15s','Location','NorthEast');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'9_3D','-dpdf','-r0')
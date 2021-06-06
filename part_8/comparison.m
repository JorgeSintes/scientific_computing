%% Comparison Van der Pol mu=1.5
clc
clear
close all

tspan = [0 40];
t0 = tspan(1);
tf = tspan(end);
h0 = 0.1;
x0 = [1;1];
mu = 1.5;
args = {mu};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [1e-2 1e-4];
slashed_tol = [2];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'ESDIRK23, tol = ' + string(tol), 'ode45, tol = ' + string(tol)];
end

method = 'ESDIRK23';
fun  = 'VanderPolFun';
jac  = 'VanderPolJac';

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    [T,X,Gout,info,stats] = ESDIRK(fun,jac,t0,tf,x0,h0,abstol,reltol,method,mu);
    X = X';
    infos(1,count) = info.nFun;
    infos(2,count) = info.nStep;
    infos(3,count) = info.nAccept;
    infos(4,count) = info.nFail;

    opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
    sol_3D_45 = ode45(@(t,x) VanderPolFun(t,x,mu),tspan, x0, opts);
    T_45 = sol_3D_45.x;
    X_45 = sol_3D_45.y;
    info_45(1) = sol_3D_45.stats.nfevals;
    info_45(3) = sol_3D_45.stats.nsteps;
    info_45(4) = sol_3D_45.stats.nfailed;
    info_45(2) = info_45(3) + info_45(4);
    infos_45(:,count) = info_45;

    subplot(2,2,1);
    if any(slashed_tol == count)
        plot(T,X(1,:),'--');
        hold on
    else
        plot(T,X(1,:));
        hold on
    end
    if any(slashed_tol == count)
        plot(T_45,X_45(1,:),'--');
        hold on
    else
        plot(T_45,X_45(1,:));
        hold on
    end
    xlabel('t');
    ylabel('x_1');

    subplot(2,2,3);
    if any(slashed_tol == count)
        plot(T,X(2,:),'--');
        hold on
    else
        plot(T,X(2,:));
        hold on
    end
    if any(slashed_tol == count)
        plot(T_45,X_45(2,:),'--');
        hold on
    else
        plot(T_45,X_45(2,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');

    subplot(2,2,[2,4]);
    if any(slashed_tol == count)
        plot(X(1,:), X(2,:),'--');
        hold on
    else
        plot(X(1,:), X(2,:));
        hold on
    end

    if any(slashed_tol == count)
        plot(X_45(1,:), X_45(2,:),'--');
        hold on
    else
        plot(X_45(1,:), X_45(2,:));
        hold on
    end

    xlabel('x_1');
    ylabel('x_2');
    xlim([-2.5 2.5])
    legend(tols_legend,'Location','SouthEast');
end
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'8_6_mu_1_5','-dpdf','-r0')

infos_csv = [tols;infos];
infos_45_csv = [tols;infos_45];
csvwrite('8_6_mu_1_5.csv',infos_csv,0,1);
csvwrite('8_6_mu_1_5_ode45.csv',infos_45_csv,0,1);

%% Comparison Van der Pol mu=15
clc
clear
close all

tspan = [0 100];
t0 = tspan(1);
tf = tspan(end);
h0 = 0.1;
x0 = [1;1];
mu = 15;
args = {mu};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [1e-2 1e-4];
slashed_tol = [2];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'ESDIRK23, tol = ' + string(tol), 'ode15s, tol = ' + string(tol)];
end

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

method = 'ESDIRK23';
fun  = 'VanderPolFun';
jac  = 'VanderPolJac';
count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    [T,X,Gout,info,stats] = ESDIRK(fun,jac,t0,tf,x0,h0,abstol,reltol,method,mu);
    X = X';
    infos(1,count) = info.nFun;
    infos(2,count) = info.nStep;
    infos(3,count) = info.nAccept;
    infos(4,count) = info.nFail;

    opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
    sol_15 = ode15s(@(t,x) VanderPolFun(t,x,mu),tspan, x0, opts);
    T_15 = sol_15.x;
    X_15 = sol_15.y;
    info_15(1) = sol_15.stats.nfevals;
    info_15(3) = sol_15.stats.nsteps;
    info_15(4) = sol_15.stats.nfailed;
    info_15(2) = info_15(3) + info_15(4);
    infos_15(:,count) = info_15;

    subplot(2,2,1);
    if any(slashed_tol == count)
        plot(T,X(1,:),'--');
        hold on
    else
        plot(T,X(1,:));
        hold on
    end
    if any(slashed_tol == count)
        plot(T_15,X_15(1,:),'--');
        hold on
    else
        plot(T_15,X_15(1,:));
        hold on
    end
    xlabel('t');
    ylabel('x_1');

    subplot(2,2,3);
    if any(slashed_tol == count)
        plot(T,X(2,:),'--');
        hold on
    else
        plot(T,X(2,:));
        hold on
    end
    if any(slashed_tol == count)
        plot(T_15,X_15(2,:),'--');
        hold on
    else
        plot(T_15,X_15(2,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');

    subplot(2,2,[2,4]);
    if any(slashed_tol == count)
        plot(X(1,:), X(2,:),'--');
        hold on
    else
        plot(X(1,:), X(2,:));
        hold on
    end

    if any(slashed_tol == count)
        plot(X_15(1,:), X_15(2,:),'--');
        hold on
    else
        plot(X_15(1,:), X_15(2,:));
        hold on
    end

    xlabel('x_1');
    ylabel('x_2');
    xlim([-2.5 2.5])
    legend(tols_legend,'Location','SouthEast');
end
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'8_6_mu_15','-dpdf','-r0')

infos_csv = [tols;infos];
infos_15_csv = [tols;infos_15];
csvwrite('8_6_mu_15.csv',infos_csv,0,1);
csvwrite('8_6_mu_15_ode15s.csv',infos_15_csv,0,1);

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
tols =  [1e-3];
slashed_tol = [3];
lim1 = [0 100];
lim2 = [0 100];
lim3 = [0 1.2];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'ESDIRK23, tol = ' + string(tol), 'ode45, tol = ' + string(tol), 'ode15s, tol = ' + string(tol)];
end

method = 'ESDIRK23';

count = 0;
infos_3D = zeros(4,size(tols,2));
infos_3D_45 = zeros(4,size(tols,2));
infos_1D = zeros(4,size(tols,2));

% figure()
% figure()
figure('Renderer', 'painters', 'Position', [680   500   900   350])
figure('Renderer', 'painters', 'Position', [680   500   900   350])
for tol=tols
    count = count+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions
    x0 = [[0; 0; Tin]];
    x0_45 = x0;
    x0_15 = x0;
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_3D = [];
    X_3D = [];
    r_3D = [];
    h_3D = [];
    info_3D = zeros(4,1);
    T_3D_45 = [];
    X_3D_45 = [];
    info_45 = zeros(4,1);
    T_3D_15 = [];
    X_3D_15 = [];
    info_15 = zeros(4,1);
    
    opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,e,struct_info_local,stats_local] = ESDIRK(@CSTR_3DFun,@CSTR_3DJac,tspan(1),tspan(2),x0,h0,abstol,reltol,method,args);
        X_local = X_local';
        T_local = T_local';
        
        info_local(1) = struct_info_local.nFun;
        info_local(2) = struct_info_local.nStep;
        info_local(3) = struct_info_local.nAccept;
        info_local(4) = struct_info_local.nFail;
        
        r_local = stats_local.r;
        h_local = stats_local.h;
        
        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];
        r_3D = [r_3D, r_local(:,1:end)];
        h_3D = [h_3D, h_local(:,1:end)];
        info_3D = info_3D + info_local';
        
        opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
        sol_3D_45 = ode45(@(t,x) CSTR_3D(t,x,F,params),tspan, x0_45, opts);
        T_3D_45 = [T_3D_45, sol_3D_45.x(1:end-1)];
        X_3D_45 = [X_3D_45, sol_3D_45.y(:,1:end-1)];
        info_45(1) = info_45(1) + sol_3D_45.stats.nfevals;
        info_45(3) = info_45(3) + sol_3D_45.stats.nsteps;
        info_45(4) = info_45(4) + sol_3D_45.stats.nfailed;
        info_45(2) = info_45(2) + info_45(3) + info_45(4);
        
        sol_3D_15 = ode15s(@(t,x) CSTR_3D(t,x,F,params),tspan, x0_15, opts);
        T_3D_15 = [T_3D_15, sol_3D_15.x(1:end-1)];
        X_3D_15 = [X_3D_15, sol_3D_15.y(:,1:end-1)];
        info_15(1) = info_15(1) + sol_3D_15.stats.nfevals;
        info_15(3) = info_15(3) + sol_3D_15.stats.nsteps;
        info_15(4) = info_15(4) + sol_3D_15.stats.nfailed;
        info_15(2) = info_15(2) + info_15(3) + info_15(4);
        
        x0 = X_local(:,end);
        x0_45 = sol_3D_45.y(:,end);
        x0_15 = sol_3D_15.y(:,end);
        h0 = h_local(end);
    end
    
    % Add last point and normalize
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    T_3D_45 = [T_3D_45,sol_3D_45.x(end)]/60;
    X_3D_45 = [X_3D_45,sol_3D_45.y(:,end)]-273;
    T_3D_15 = [T_3D_15,sol_3D_15.x(end)]/60;
    X_3D_15 = [X_3D_15,sol_3D_15.y(:,end)]-273;
    r_3D = [r_3D,r_local(end)];
    h_3D = [h_3D, h_local(:,end)]/60;
    infos_3D(:,count) = info_3D;
    infos_3D_45(:,count) = info_45;
    infos_3D_15(:,count) = info_15;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions
    x0 = Tin;
    x0_45 = x0;
    x0_15 = x0;
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_1D = [];
    X_1D = [];
    r_1D = [];
    h_1D = [];
    info_1D = zeros(4,1);
    T_1D_45 = [];
    X_1D_45 = [];
    info_45 = zeros(4,1);
    T_1D_15 = [];
    X_1D_15 = [];
    info_15 = zeros(4,1);
    
    opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        params = [beta,k0,EaR,CAin,CBin,Tin,V];
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,e,struct_info_local,stats_local] = ESDIRK(@CSTR_1DFun,@CSTR_1DJac,tspan(1),tspan(2),x0,h0,abstol,reltol,method,args);
        
        X_local = X_local';
        T_local = T_local';
        
        info_local(1) = struct_info_local.nFun;
        info_local(2) = struct_info_local.nStep;
        info_local(3) = struct_info_local.nAccept;
        info_local(4) = struct_info_local.nFail;
        
        r_local = stats_local.r;
        h_local = stats_local.h;
        
        T_1D = [T_1D, T_local(1:end-1)];
        X_1D = [X_1D, X_local(:,1:end-1)];
        r_1D = [r_1D, r_local(:,1:end)];
        h_1D = [h_1D, h_local(:,1:end)];
        info_1D = info_1D + info_local';
        
        opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
        sol_1D_45 = ode45(@(t,x) CSTR_1D(t,x,F,params),tspan, x0_45, opts);
        T_1D_45 = [T_1D_45, sol_1D_45.x(1:end-1)];
        X_1D_45 = [X_1D_45, sol_1D_45.y(:,1:end-1)];
        info_45(1) = info_45(1) + sol_1D_45.stats.nfevals;
        info_45(3) = info_45(3) + sol_1D_45.stats.nsteps;
        info_45(4) = info_45(4) + sol_1D_45.stats.nfailed;
        info_45(2) = info_45(2) + info_45(3) + info_45(4);
        
        sol_1D_15 = ode15s(@(t,x) CSTR_1D(t,x,F,params),tspan, x0_15, opts);
        T_1D_15 = [T_1D_15, sol_1D_15.x(1:end-1)];
        X_1D_15 = [X_1D_15, sol_1D_15.y(:,1:end-1)];
        info_15(1) = info_15(1) + sol_1D_15.stats.nfevals;
        info_15(3) = info_15(3) + sol_1D_15.stats.nsteps;
        info_15(4) = info_15(4) + sol_1D_15.stats.nfailed;
        info_15(2) = info_15(2) + info_15(3) + info_15(4);
        
        x0 = X_local(:,end);
        x0_45 = sol_1D_15.y(:,end);
        x0_15 = sol_1D_15.y(:,end);
        h0 = h_local(end);
    end
    
    % Add last point and normalize
    T_1D = [T_1D,T_local(end)]/60;
    X_1D = [X_1D, X_local(:,end)]-273;
    T_1D_45 = [T_1D_45,sol_1D_45.x(end)]/60;
    X_1D_45 = [X_1D_45,sol_1D_45.y(:,end)]-273;
    T_1D_15 = [T_1D_15,sol_1D_15.x(end)]/60;
    X_1D_15 = [X_1D_15,sol_1D_15.y(:,end)]-273;
    r_1D = [r_1D,r_local(end)];
    h_1D = [h_1D, h_local(:,end)]/60;
    infos_1D(:,count) = info_1D;
    infos_1D_45(:,count) = info_45;
    infos_1D_15(:,count) = info_15;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTTING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(1)
    if any(slashed_tol == count)
        plot(T_3D, X_3D(3,:),'--');
        hold on
    else
        plot(T_3D, X_3D(3,:));
        hold on
    end
    if any(slashed_tol == count)
        plot(T_3D_45, X_3D_45(3,:),'--');
        hold on
    else
        plot(T_3D_45, X_3D_45(3,:));
        hold on
    end
    if any(slashed_tol == count)
        plot(T_3D_15, X_3D_15(3,:),'--');
        hold on
    else
        plot(T_3D_15, X_3D_15(3,:));
        hold on
    end
    
    legend(tols_legend,'Location','NorthWest');
    ylabel('T(°C)');
    xlabel('t (min)');
    ylim(lim1)
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'8_6_3D','-dpdf','-r0')
    
    figure(2)
    if any(slashed_tol == count)
        plot(T_1D, X_1D,'--');
        hold on
    else
        plot(T_1D, X_1D);
        hold on
    end
    if any(slashed_tol == count)
        plot(T_1D_45, X_1D_45,'--');
        hold on
    else
        plot(T_1D_45, X_1D_45);
        hold on
    end
    if any(slashed_tol == count)
        plot(T_1D_15, X_1D_15,'--');
        hold on
    else
        plot(T_1D_15, X_1D_15);
        hold on
    end
    
    legend(tols_legend,'Location','NorthWest');
    ylabel('T(°C)');
    xlabel('t (min)');
    ylim(lim2)
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'8_6_1D','-dpdf','-r0')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT FILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    infos_3D_csv = [tols;infos_3D];
    infos_3D_45_csv = [tols;infos_3D_45];
    infos_3D_15_csv = [tols;infos_3D_15];
    csvwrite('8_6_3D.csv',infos_3D_csv,0,1);
    csvwrite('8_6_3D_45.csv',infos_3D_45_csv,0,1);
    csvwrite('8_6_3D_15.csv',infos_3D_15_csv,0,1);
    
    infos_1D_csv = [tols;infos_1D];
    infos_1D_45_csv = [tols;infos_1D_45];
    infos_1D_15_csv = [tols;infos_1D_15];
    csvwrite('8_6_1D.csv',infos_1D_csv,0,1);
    csvwrite('8_6_1D_45.csv',infos_1D_45_csv,0,1);
    csvwrite('8_6_1D_15.csv',infos_1D_15_csv,0,1);
end

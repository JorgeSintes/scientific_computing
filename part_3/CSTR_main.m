%% CSTR 3D Implicit Euler First Trial
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

% Initial conditions
x0 = [[0; 0; Tin]];
h = 0.5;

T = [];
X = [];

for i=1:size(Fs,2)
    tspan = tspans(i,:)*60;
    F = Fs(i);
    args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
    [T_local,X_local] = EulerImplicit_fixed(@CSTR_3D_fJ, tspan, h, x0, args);
    
    T = [T, T_local(1:end-1)];
    X = [X, X_local(:,1:end-1)];
    
    x0 = X_local(:,end);
end
% Add last point and normalize
T = [T,T_local(end)]/60;
X = [X, X_local(:,end)]-273;

T_3D = T;
X_3D = X;

figure()
subplot(2,1,1)
plot(T, X(3,:));
ylabel('T(°C)');
subplot(2,1,2)
plot(tspans_init,F_init);
ylim([100 800]);
xlabel('t (min)');
ylabel('F (mL/min)');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'3_5_3D_trial','-dpdf','-r0')


%% Flow plot
figure('Renderer', 'painters', 'Position', [10 10 900 300])
plot(tspans_init,F_init);
ylim([100 800]);
xlabel('t (min)');
ylabel('F (mL/min)');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'3_5_flow','-dpdf','-r0')

%% CSTR 1D Explicit Euler First Trial
clc
% clearvars -except T_3D X_3D
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
    F_init = [F_init,Fs(i),Fs(i)];   
end

% Initial conditions
x0 = Tin;
h = 0.5;
tspan = [0 35];

T = [];
X = [];

for i=1:size(Fs,2)
    tspan = tspans(i,:)*60;
    F = Fs(i);
    args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
    [T_local,X_local] = EulerImplicit_fixed(@CSTR_1D_fJ, tspan, h, x0, args);
    
    T = [T, T_local(1:end-1)];
    X = [X, X_local(:,1:end-1)];
    
    x0 = X_local(:,end);
end
% Add last point and normalize
T = [T,T_local(end)]/60;
X = [X, X_local(:,end)]-273;

figure()
subplot(2,1,1)
plot(T, X);
ylabel('T(°C)');

subplot(2,1,2)
plot(tspans_init,F_init);
ylim([100 800]);
xlabel('t (min)');
ylabel('F (mL/min)');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'3_5_1D_trial','-dpdf','-r0')

%% Plot 3D against 1D
figure()
subplot(2,1,1)
plot(T_3D, X_3D(3,:));
hold on
plot(T, X, '--');
ylabel('T(°C)');
legend('3D', '1D');

subplot(2,1,2)
plot(T_3D, X_3D(3,:));
hold on
plot(T, X, '--');
xlabel('t (min)');
ylabel('T(°C)');
legend('3D', '1D');
xlim([0 1.5]);
ylim([0 5]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'3_5_3D_vs_1D','-dpdf','-r0')

%% CSTR 3D-1D explicit fixed for different timesteps
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
% TIME-STEP SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hs = [0.15,0.1,0.05,0.01];
slashed_h = [4];

h_legend = [];
for h=hs
    h_legend = [h_legend, 'h = '+string(h)];
end

count = 0;

% figure('Renderer', 'painters', 'Position', [10 10 1200 600])
figure('Renderer', 'painters', 'Position', [680   678   560   210])
figure('Renderer', 'painters', 'Position', [680   678   560   210])
for h=hs*60
    count = count+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions
    x0 = [[0; 0; Tin]];
    tspan = [0 35];

    T_3D = [];
    X_3D = [];

    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local] = EulerImplicit_fixed(@CSTR_3D_fJ, tspan, h, x0, args);

        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];

        x0 = X_local(:,end);
    end
    % Add last point and normalize
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions 1D
    x0 = Tin;
    tspan = [0 35];

    T_1D = [];
    X_1D = [];

    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local] = EulerImplicit_fixed(@CSTR_1D_fJ, tspan, h, x0, args);

        T_1D = [T_1D, T_local(1:end-1)];
        X_1D = [X_1D, X_local(:,1:end-1)];

        x0 = X_local(:,end);
    end
    % Add last point and normalize
    T_1D = [T_1D,T_local(end)]/60;
    X_1D = [X_1D, X_local(:,end)]-273;
    
    figure(1)
    if any(slashed_h == count)
        plot(T_3D, X_3D(3,:),'--');
        hold on
    else
        plot(T_3D, X_3D(3,:));
        hold on
    end
    
    legend(h_legend,'Location','NorthWest');
    ylabel('T(°C)');
    xlabel('t (min)');
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'3_5_3D_hs','-dpdf','-r0')
    
    figure(2)
    if any(slashed_h == count)
        plot(T_1D, X_1D,'--');
        hold on
    else
        plot(T_1D, X_1D);
        hold on
    end
    ylabel('T(°C)');
    xlabel('t (min)');
    legend(h_legend,'Location','NorthWest');
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'3_5_1D_hs','-dpdf','-r0')

end


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
tols =  [2e-2 5e-3 1e-3 1e-5];
slashed_tol = [4];
lim1 = [0 120];
lim2 = [-0.02 2.5];
lim3 = [0 1.2];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'tol = ' + string(tol)];
end

count = 0;
infos_3D = zeros(4,size(tols,2));
infos_1D = zeros(4,size(tols,2));

% figure()
% figure()
figure('Renderer', 'painters', 'Position', [680   500   900   600])
figure('Renderer', 'painters', 'Position', [680   500   900   600])
for tol=tols
    count = count+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions
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
    % Add last point and normalize
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    r_3D = [r_3D,r_local(end)];
    h_3D = [h_3D, h_local(:,end)]/60;
    infos_3D(:,count) = info_3D;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions
    x0 = Tin;
    tspan = [0 35];
    h0 = 0.01*60;
    abstol = tol;
    reltol = tol;
    
    T_1D = [];
    X_1D = [];
    r_1D = [];
    h_1D = [];
    info_1D = zeros(4,1);
    
    for i=1:size(Fs,2)
        tspan = tspans(i,:)*60;
        F = Fs(i);
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
        [T_local,X_local,r_local,h_local, info_local] = EulerImplicit_adaptive(@CSTR_1D_fJ,tspan,h0,x0,abstol,reltol,args);

        T_1D = [T_1D, T_local(1:end-1)];
        X_1D = [X_1D, X_local(:,1:end-1)];
        r_1D = [r_1D, r_local(:,1:end)];
        h_1D = [h_1D, h_local(:,1:end)];
        info_1D = info_1D + info_local';
        
        x0 = X_local(:,end);
        h0 = h_local(end);
    end
    % Add last point and normalize
    T_1D = [T_1D,T_local(end)]/60;
    X_1D = [X_1D, X_local(:,end)]-273;
    r_1D = [r_1D,r_local(end)];
    h_1D = [h_1D, h_local(:,end)]/60;
    infos_1D(:,count) = info_1D;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTTING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(1)
    subplot(3,1,1)
    if any(slashed_tol == count)
        plot(T_3D, X_3D(3,:),'--');
        hold on
    else
        plot(T_3D, X_3D(3,:));
        hold on
    end
    legend(tols_legend,'Location','NorthWest');
    ylabel('T(°C)');
    xlabel('t (min)');
    ylim(lim1)
    
    subplot(3,1,2)
    if any(slashed_tol == count)
        plot(T_3D, h_3D,'--');
        hold on
    else
        plot(T_3D, h_3D);
        hold on
    end
    xlabel('t (min)');
    ylabel('h (min)');
    ylim(lim2)
    
    subplot(3,1,3)
    if any(slashed_tol == count)
        plot(T_3D, r_3D,'--');
        hold on
    else
        plot(T_3D, r_3D);
        hold on
    end
    plot([0 35], [1 1],'r');
    xlabel('t (min)');
    ylabel('r');
    ylim(lim3)
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'3_5_3D_tols','-dpdf','-r0')
  
    figure(2)
    subplot(3,1,1)
    if any(slashed_tol == count)
        plot(T_1D, X_1D,'--');
        hold on
    else
        plot(T_1D, X_1D);
        hold on
    end
    legend(tols_legend,'Location','NorthWest');
    ylabel('T(°C)');
    xlabel('t (min)');
    ylim(lim1)
    
    subplot(3,1,2)
    if any(slashed_tol == count)
        plot(T_1D, h_1D,'--');
        hold on
    else
        plot(T_1D, h_1D);
        hold on
    end
    xlabel('t (min)');
    ylabel('h (min)');
    ylim(lim2)
    
    subplot(3,1,3)
    if any(slashed_tol == count)
        plot(T_1D, r_1D,'--');
        hold on
    else
        plot(T_1D, r_1D);
        hold on
    end
    plot([0 35], [1 1],'r');
    xlabel('t (min)');
    ylabel('r');
    ylim(lim3)
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'3_5_1D_tols','-dpdf','-r0')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT FILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    infos_3D_csv = [tols;infos_3D];
    csvwrite('3_5_CSTR_3D.csv',infos_3D_csv,0,1);
    
    infos_1D_csv = [tols;infos_1D];
    csvwrite('3_5_CSTR_1D.csv',infos_1D_csv,0,1);
    

end

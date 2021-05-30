%% CSTR 3D Explicit Euler First Trial
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
tspan = [0 35];

T = [];
X = [];

for i=1:size(Fs,2)
    tspan = tspans(i,:)*60;
    F = Fs(i);
    args = {F, [beta,k0,EaR,CAin,CBin,Tin,V]};
    [T_local,X_local] = EulerExplicit_fixed(@CSTR_3D, tspan, h, x0, args);
    
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
plot(T, X(3,:)-273);
ylabel('T(°C)');
subplot(2,1,2)
plot(tspans_init,F_init);
ylim([100 800]);
xlabel('t (min)');
ylabel('F (mL/min)');


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
    F_init = [F_init,Fs(i), Fs(i)];   
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
    [T_local,X_local] = EulerExplicit_fixed(@CSTR_1D, tspan, h, x0, args);
    
    T = [T, T_local(1:end-1)];
    X = [X, X_local(:,1:end-1)];
    
    x0 = X_local(:,end);
end
% Add last point and normalize
T = [T,T_local(end)]/60;
X = [X, X_local(:,end)]-273;

figure()
subplot(2,1,1)
plot(T_3D, X_3D(3,:));
hold on
plot(T, X);
ylabel('T(°C)');
legend('3D', '1D');

subplot(2,1,2)
plot(tspans_init,F_init);
ylim([100 800]);
xlabel('t (min)');
ylabel('F (mL/min)');


%% Plot one against other


%%
for h=hs
    [T,X] = EulerExplicit_fixed(@VanderPol, tspan, h, x0, args);

    figure()
    plot(T,X(1,:));
    hold on
    plot(T,X(2,:));
    xlabel('time')
    ylabel('x(t)')
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'2_4_vdp_mu_1_5_h_' + strrep(string(h),'.','_')+'_T','-dpdf','-r0')
    
    figure()
    plot(X(1,:), X(2,:));
    xlabel('x_1');
    ylabel('x_2');
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'2_4_fixed_vdp_mu_1_5_h_' + strrep(string(h),'.','_')+'_PP','-dpdf','-r0')
end

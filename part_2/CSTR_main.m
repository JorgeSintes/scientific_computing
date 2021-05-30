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
print(gcf,'2_5_CSTR_3D_trial','-dpdf','-r0')


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
print(gcf,'2_5_CSTR_1D_trial','-dpdf','-r0')

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
print(gcf,'2_5_CSTR_3D_vs_1D','-dpdf','-r0')

%% CSTR 3D-1D explicit for different timesteps
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
slashed_h = [3,4];

h_legend = [];
for h=hs
    h_legend = [h_legend, 'h = '+string(h)];
end

count = 0;

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
        [T_local,X_local] = EulerExplicit_fixed(@CSTR_3D, tspan, h, x0, args);

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
        [T_local,X_local] = EulerExplicit_fixed(@CSTR_1D, tspan, h, x0, args);

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
    legend(h_legend);
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'2_5_CSTR_3D_hs','-dpdf','-r0')
    
    figure(2)
    if any(slashed_h == count)
        plot(T_1D, X_1D,'--');
        hold on
    else
        plot(T_1D, X_1D);
        hold on
    end
    ylabel('T(°C)');
    legend(h_legend);
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'2_5_CSTR_1D_hs','-dpdf','-r0')

end

%% CSTR 3D Euler Maruyama
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
h = 0.1;

Ns = 5;
sigma_EM_3D = 5;
sigma_IE_3D = 5;
sigma_EM_1D = 5;
sigma_IE_1D = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
x0 = [[0; 0; Tin]];
tspan = [0 35];

X = [];

for i=1:Ns
    T_3D = [];
    X_3D = [];
    for j=1:size(Fs,2)
        tspan = tspans(j,:)*60;
        F = Fs(j);
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], sigma_EM_3D};
        [T_local,X_local] = EulerMaruyama(@CSTR_3D_Drift, @CSTR_3D_Diffusion, tspan, h, x0, 1, args);

        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];

        x0 = X_local(:,end);
    end
    % Add last point and normalize
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    X(:,:,i) = X_3D;
end

T_3D = [];
X_3D = [];
for i=1:size(Fs,2)
    tspan = tspans(i,:)*60;
    F = Fs(i);
    args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], 0};
    [T_local,X_local] = EulerMaruyama(@CSTR_3D_Drift, @CSTR_3D_Diffusion, tspan, h, x0, 1, args);

    T_3D = [T_3D, T_local(1:end-1)];
    X_3D = [X_3D, X_local(:,1:end-1)];

    x0 = X_local(:,end);
end
% Add last point and normalize
T_3D = [T_3D,T_local(end)]/60;
X_3D = [X_3D, X_local(:,end)]-273;

Xd = X_3D;

figure('Renderer', 'painters', 'Position', [680   678   800   300])

for i=1:Ns
    plot(T_3D, X(3,:,i));
    hold on
end
plot(T_3D, Xd(3,:),'LineWidth',1.5,'Color','black');
ylabel('T(°C)');
xlabel('t (min)');
ylim([0 90]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_6_EulerMaruyama_3D','-dpdf','-r0')

%% CSTR 3D Implicit-Explicit
clc
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
x0 = [[0; 0; Tin]];
tspan = [0 35];

X = [];

for i=1:Ns
    T_3D = [];
    X_3D = [];
    for j=1:size(Fs,2)
        tspan = tspans(j,:)*60;
        F = Fs(j);
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], sigma_IE_3D};
        [T_local,X_local] = ImplicitExplicit(@CSTR_3D_Drift, @CSTR_3D_Diffusion, tspan, h, x0, 1, args);

        T_3D = [T_3D, T_local(1:end-1)];
        X_3D = [X_3D, X_local(:,1:end-1)];

        x0 = X_local(:,end);
    end
    % Add last point and normalize
    T_3D = [T_3D,T_local(end)]/60;
    X_3D = [X_3D, X_local(:,end)]-273;
    
    X(:,:,i) = X_3D;
end

T_3D = [];
X_3D = [];
for i=1:size(Fs,2)
    tspan = tspans(i,:)*60;
    F = Fs(i);
    args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], 0};
    [T_local,X_local] = ImplicitExplicit(@CSTR_3D_Drift, @CSTR_3D_Diffusion, tspan, h, x0, 1, args);

    T_3D = [T_3D, T_local(1:end-1)];
    X_3D = [X_3D, X_local(:,1:end-1)];

    x0 = X_local(:,end);
end
% Add last point and normalize
T_3D = [T_3D,T_local(end)]/60;
X_3D = [X_3D, X_local(:,end)]-273;

Xd = X_3D;

figure('Renderer', 'painters', 'Position', [680   678   800   300])

for i=1:Ns
    plot(T_3D, X(3,:,i));
    hold on
end
plot(T_3D, Xd(3,:),'LineWidth',1.5,'Color','black');
ylabel('T(°C)');
xlabel('t (min)');
ylim([0 90]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_6_ImplicitExplicit_3D','-dpdf','-r0')

%% CSTR 1D Euler Maruyama
clc
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
x0 = Tin;
tspan = [0 35];

X = [];

for i=1:Ns
    T_1D = [];
    X_1D = [];
    for j=1:size(Fs,2)
        tspan = tspans(j,:)*60;
        F = Fs(j);
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], sigma_EM_1D};
        [T_local,X_local] = EulerMaruyama(@CSTR_1D_Drift, @CSTR_1D_Diffusion, tspan, h, x0, 1, args);

        T_1D = [T_1D, T_local(1:end-1)];
        X_1D = [X_1D, X_local(:,1:end-1)];

        x0 = X_local(:,end);
    end
    % Add last point and normalize
    T_1D = [T_1D,T_local(end)]/60;
    X_1D = [X_1D, X_local(:,end)]-273;
    
    X(:,i) = X_1D;
end

T_1D = [];
X_1D = [];
for i=1:size(Fs,2)
    tspan = tspans(i,:)*60;
    F = Fs(i);
    args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], 0};
    [T_local,X_local] = EulerMaruyama(@CSTR_1D_Drift, @CSTR_1D_Diffusion, tspan, h, x0, 1, args);

    T_1D = [T_1D, T_local(1:end-1)];
    X_1D = [X_1D, X_local(:,1:end-1)];

    x0 = X_local(:,end);
end
% Add last point and normalize
T_1D = [T_1D,T_local(end)]/60;
X_1D = [X_1D, X_local(:,end)]-273;

Xd = X_1D;

figure('Renderer', 'painters', 'Position', [680   678   800   300])

for i=1:Ns
    plot(T_1D, X(:,i));
    hold on
end
plot(T_1D, Xd,'LineWidth',1.5,'Color','black');
ylabel('T(°C)');
xlabel('t (min)');
ylim([0 90]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_6_EulerMaruyama_1D','-dpdf','-r0')

%% CSTR 1D Euler Maruyama
clc
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
x0 = Tin;
tspan = [0 35];

X = [];

for i=1:Ns
    T_1D = [];
    X_1D = [];
    for j=1:size(Fs,2)
        tspan = tspans(j,:)*60;
        F = Fs(j);
        args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], sigma_IE_1D};
        [T_local,X_local] = ImplicitExplicit(@CSTR_1D_Drift, @CSTR_1D_Diffusion, tspan, h, x0, 1, args);

        T_1D = [T_1D, T_local(1:end-1)];
        X_1D = [X_1D, X_local(:,1:end-1)];

        x0 = X_local(:,end);
    end
    % Add last point and normalize
    T_1D = [T_1D,T_local(end)]/60;
    X_1D = [X_1D, X_local(:,end)]-273;
    
    X(:,i) = X_1D;
end

T_1D = [];
X_1D = [];
for i=1:size(Fs,2)
    tspan = tspans(i,:)*60;
    F = Fs(i);
    args = {F, [beta,k0,EaR,CAin,CBin,Tin,V], 0};
    [T_local,X_local] = ImplicitExplicit(@CSTR_1D_Drift, @CSTR_1D_Diffusion, tspan, h, x0, 1, args);

    T_1D = [T_1D, T_local(1:end-1)];
    X_1D = [X_1D, X_local(:,1:end-1)];

    x0 = X_local(:,end);
end
% Add last point and normalize
T_1D = [T_1D,T_local(end)]/60;
X_1D = [X_1D, X_local(:,end)]-273;

Xd = X_1D;

figure('Renderer', 'painters', 'Position', [680   678   800   300])

for i=1:Ns
    plot(T_1D, X(:,i));
    hold on
end
plot(T_1D, Xd,'LineWidth',1.5,'Color','black');
ylabel('T(°C)');
xlabel('t (min)');
ylim([0 90]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_6_ImplicitExplicit_1D','-dpdf','-r0')


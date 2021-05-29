%% 2.2 Explicit Fixed Van der Pol mu=1.5
clc
clear
close all

x0 = [1;1];
hs = [0.2 0.1 0.01 0.001];

tspan = [0 40];
mu = 1.5;
args = {mu};

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

%% 2.2 Explicit Fixed Van der Pol mu=15
clc
clear
close all

x0 = [1;1];
hs = [0.05 0.01 0.001];

tspan = [0 100];
mu = 15;
args = {mu};

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
    print(gcf,'2_4_fixed_vdp_mu_15_h_' + strrep(string(h),'.','_')+'_T','-dpdf','-r0')
    
    figure()
    plot(X(1,:), X(2,:));
    xlabel('x_1');
    ylabel('x_2');
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'2_4_vdp_mu_15_h_' + strrep(string(h),'.','_')+'_PP','-dpdf','-r0')
end

%% 2.3. Explicit Adaptive Van der Pol
clc
clear
close all




%% Explicit Fixed Van der Pol mu=1.5
clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME-STEP SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hs = [0.2 0.1 0.01 0.001];
slashed_h = [4];


h_legend = [];
for h=hs
    h_legend = [h_legend, 'h = '+string(h)];
end

x0 = [1;1];
tspan = [0 40];
mu = 1.5;
args = {mu};

figure('Renderer', 'painters', 'Position', [10 10 1200 533])

count = 0;
for h=hs
    count = count+1;
    [T,X] = EulerExplicit_fixed(@VanderPol, tspan, h, x0, args);
    
    subplot(2,2,1)
    if any(slashed_h == count)
        plot(T,X(1,:),'--');
        hold on
    else
        plot(T,X(1,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,3)
    if any(slashed_h == count)
        plot(T,X(2,:),'--');
        hold on
    else
        plot(T,X(2,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,[2,4])
    if any(slashed_h == count)
        plot(X(1,:), X(2,:),'--');
        hold on
    else
        plot(X(1,:), X(2,:));
        hold on
    end
    xlabel('x_1');
    ylabel('x_2');
    legend(h_legend,'Location','SouthEast')
end

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'2_4_fixed_mu_1_5','-dpdf','-r0')

%% Explicit Fixed Van der Pol mu=15
clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME-STEP SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hs = [0.01 0.005 0.001 0.0001];
slashed_h = [4];


h_legend = [];
for h=hs
    h_legend = [h_legend, 'h = '+string(h)];
end

x0 = [1;1];
tspan = [0 100];
mu = 15;
args = {mu};

figure('Renderer', 'painters', 'Position', [10 10 1200 533])

count = 0;
for h=hs
    count = count+1;
    [T,X] = EulerExplicit_fixed(@VanderPol, tspan, h, x0, args);
    
    subplot(2,2,1)
    if any(slashed_h == count)
        plot(T,X(1,:),'--');
        hold on
    else
        plot(T,X(1,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,3)
    if any(slashed_h == count)
        plot(T,X(2,:),'--');
        hold on
    else
        plot(T,X(2,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,[2,4])
    if any(slashed_h == count)
        plot(X(1,:), X(2,:),'--');
        hold on
    else
        plot(X(1,:), X(2,:));
        hold on
    end
    xlabel('x_1');
    ylabel('x_2');
    legend(h_legend,'Location','SouthEast')
end

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'2_4_fixed_mu_15','-dpdf','-r0')

%% Explicit Adaptive Van der Pol mu=1.5
clc
clear
close all

tspan = [0 40];
h0 = 0.1;
x0 = [1;1];
mu = 1.5;
args = {mu};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [5e-2, 5e-3, 1e-4, 1e-6];
slashed_tol = [4];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'tol = ' + string(tol)];
end

figure('Renderer', 'painters', 'Position', [10 10 1200 800])

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;
    
    [T,X,r_out,h_out,info] = EulerExplicit_adaptive(@VanderPol, tspan, h0, x0, abstol, reltol, args);
    infos(:,count) = info;

    subplot(3,2,1);
    if any(slashed_tol == count)
        plot(T,X(1,:),'--');
        hold on
    else
        plot(T,X(1,:));
        hold on
    end
    xlabel('t');
    ylabel('x_1');

    subplot(3,2,3);
    if any(slashed_tol == count)
        plot(T,X(2,:),'--');
        hold on
    else
        plot(T,X(2,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');

    subplot(3,2,[2,4]);
    if any(slashed_tol == count)
        plot(X(1,:), X(2,:),'--');
        hold on
    else
        plot(X(1,:), X(2,:));
        hold on
    end
    xlabel('x_1');
    ylabel('x_2');
    xlim([-2.5 2.5])
    legend([tols_legend],'Location','SouthEast');

    subplot(3,2,5);
    if any(slashed_tol == count)
        plot(T(1:end-1),h_out,'--');
        hold on
    else
        plot(T(1:end-1),h_out);
        hold on
    end
    xlabel('t');
    ylabel('h');
    ylim([-0.02 0.35])

    subplot(3,2,6);
    if any(slashed_tol == count)
        plot(T(1:end-1),r_out,'--');
        hold on
    else
        plot(T(1:end-1),r_out);
        hold on
    end
    plot(tspan, [1 1],'r');
    xlabel('t');
    ylabel('r');
    ylim([0 1.2]);
end

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'2_4_adaptive_mu_1_5','-dpdf','-r0')

infos_csv = [tols;infos];
csvwrite('2_4_adaptive_mu_1_5.csv',infos_csv,0,1);

%% Explicit Adaptive Van der Pol mu=15
clc
clear
close all

tspan = [0 100];
h0 = 0.1;
x0 = [1;1];
mu = 15;
args = {mu};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [5e-2, 5e-3, 1e-4, 1e-6];
slashed_tol = [4];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'tol = ' + string(tol)];
end

figure('Renderer', 'painters', 'Position', [10 10 1200 800])

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;
    
    [T,X,r_out,h_out,info] = EulerExplicit_adaptive(@VanderPol, tspan, h0, x0, abstol, reltol, args);
    infos(:,count) = info;

    subplot(3,2,1);
    if any(slashed_tol == count)
        plot(T,X(1,:),'--');
        hold on
    else
        plot(T,X(1,:));
        hold on
    end
    xlabel('t');
    ylabel('x_1');

    subplot(3,2,3);
    if any(slashed_tol == count)
        plot(T,X(2,:),'--');
        hold on
    else
        plot(T,X(2,:));
        hold on
    end
    xlabel('t');
    ylabel('x_2');

    subplot(3,2,[2,4]);
    if any(slashed_tol == count)
        plot(X(1,:), X(2,:),'--');
        hold on
    else
        plot(X(1,:), X(2,:));
        hold on
    end
    xlabel('x_1');
    ylabel('x_2');
    xlim([-2.5 2.5])
    legend([tols_legend],'Location','SouthEast');

    subplot(3,2,5);
    if any(slashed_tol == count)
        plot(T(1:end-1),h_out,'--');
        hold on
    else
        plot(T(1:end-1),h_out);
        hold on
    end
    xlabel('t');
    ylabel('h');
    ylim([-0.02 0.25])

    subplot(3,2,6);
    if any(slashed_tol == count)
        plot(T(1:end-1),r_out,'--');
        hold on
    else
        plot(T(1:end-1),r_out);
        hold on
    end
    plot(tspan, [1 1],'r');
    xlabel('t');
    ylabel('r');
    ylim([0 1.2]);
end


set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'2_4_adaptive_mu_15','-dpdf','-r0')

infos_csv = [tols;infos];
csvwrite('2_4_adaptive_mu_15.csv',infos_csv,0,1);


%% Explicit Fixed Van der Pol mu=1.5
clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME-STEP SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hs = [0.5 0.25 0.1 0.01];
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
    [T,X] = ClassicRK4_fixed(@VanderPol, tspan, h, x0, args);
    
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
print(gcf,'5_4_RK4_mu_1_5','-dpdf','-r0')

%% Explicit Fixed Van der Pol mu=15
clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME-STEP SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hs = [0.05 0.025 0.01 0.001];
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
    [T,X] = ClassicRK4_fixed(@VanderPol, tspan, h, x0, args);
    
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
print(gcf,'5_4_RK4_mu_15','-dpdf','-r0')


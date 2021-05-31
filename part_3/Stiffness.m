%% Explicit Adaptive Van der Pol mu=15
clc
clear
close all

tspan = [0 200];
h0 = 0.1;
x0 = [1;1];
mu = 100;
args = {mu};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [1e-3];
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
%     ylim([-0.02 0.25])

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
print(gcf,'3_7_explicit','-dpdf','-r0')

infos_csv = [tols;infos];
csvwrite('3_7_explicit.csv',infos_csv,0,1);



%% Implicit Adaptive Van der Pol mu=12

tspan = [0 200];
h0 = 0.1;
x0 = [0.5;0.5];
mu = 100;
args = {mu};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [1e-3];
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
    
    [T,X,r_out,h_out,info] = EulerImplicit_adaptive(@VanderPol_fJ, tspan, h0, x0, abstol, reltol, args);
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
%     ylim([-0.02 0.25])

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
print(gcf,'3_7_implicit','-dpdf','-r0')

infos_csv = [tols;infos];
csvwrite('3_7_implicit.csv',infos_csv,0,1);
%% ESDIRK23 Adaptive Van der Pol mu=1.5
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
tols =  [5e-2, 1e-2, 1e-3, 1e-5];
slashed_tol = [4];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'tol = ' + string(tol)];
end

method = 'ESDIRK23';
fun  = 'VanderPolFun';
jac  = 'VanderPolJac';

figure('Renderer', 'painters', 'Position', [10 10 1200 800])

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

    T_out = stats.t;
    h_out = stats.h;
    r_out = stats.r;
    
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
        plot(T_out,h_out,'--');
        hold on
    else
        plot(T_out,h_out);
        hold on
    end
    xlabel('t');
    ylabel('h');
    ylim([-0.02 1])

    subplot(3,2,6);
    if any(slashed_tol == count)
        plot(T_out,r_out,'--');
        hold on
    else
        plot(T_out,r_out);
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
print(gcf,'8_4_adaptive_mu_1_5','-dpdf','-r0')

infos_csv = [tols;infos];
csvwrite('8_4_adaptive_mu_1_5.csv',infos_csv,0,1);

%% ESDIRK23 Adaptive Van der Pol mu=15
clc
clear
close all

tspan = [0 40];
t0 = tspan(1);
tf = tspan(end);
h0 = 0.1;
x0 = [1;1];
mu = 15;
args = {mu};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCE SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tols =  [5e-2, 1e-2, 1e-3, 1e-5];
% tols =  [5e-2, 5e-3, 1e-4, 1e-6];

slashed_tol = [4];

tols_legend = [];

for tol=tols
    tols_legend = [tols_legend, 'tol = ' + string(tol)];
end

method = 'ESDIRK23';
fun  = 'VanderPolFun';
jac  = 'VanderPolJac';

figure('Renderer', 'painters', 'Position', [10 10 1200 800])

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

    T_out = stats.t;
    h_out = stats.h;
    r_out = stats.r;
    
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
        plot(T_out,h_out,'--');
        hold on
    else
        plot(T_out,h_out);
        hold on
    end
    xlabel('t');
    ylabel('h');
    ylim([-0.02 1])

    subplot(3,2,6);
    if any(slashed_tol == count)
        plot(T_out,r_out,'--');
        hold on
    else
        plot(T_out,r_out);
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
print(gcf,'8_4_adaptive_mu_15','-dpdf','-r0')

infos_csv = [tols;infos];
csvwrite('8_4_adaptive_mu_15.csv',infos_csv,0,1);

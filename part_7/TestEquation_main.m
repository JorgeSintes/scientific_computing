%% 1.3.Test Equation error measurements
clc
clear
close all

tspan = [0 10];
h0 = 0.2;
x0 = 1;
lambda = -1;

tol = 1e-5;

reltol = tol;
abstol = tol;

args = {lambda};

butcher = ERKSolverErrorEstimationParameters('DOPRI54');

[T,X,X_real,e_l,e_g,e_out,r_out,h_out,info] = DOPRI54TE_adaptive(@TestEquation, ...
    @Exponential, tspan, h0, x0,abstol,reltol, butcher, args);

figure('Renderer', 'painters', 'Position', [680   300   600   700])
subplot(5,1,[1,2])
yyaxis left
plot(T,X)
hold on
plot(T,X_real,'--','Color',[0.9290, 0.6940, 0.1250])

xlabel('t')
ylabel('x(t)')

yyaxis right
plot(T(1:(end-1)),e_l(1:(end-1)))
plot(T(1:(end-1)),e_g(1:(end-1)))
hold off
ylabel('errors')
legend('numerical solution','real solution', 'local error', 'global error','Location','SouthEast')
limits = ylim;

subplot(5,1,3)
plot(T(1:(end-1)),e_l(1:(end-1)))
hold on
plot(T(1:(end-1)),e_out(1:(end-1)))
ylabel('error')
legend('Real local error', 'Estimated local error')

subplot(5,1,4)
plot(T(1:(end-1)), h_out)
ylabel('h')

subplot(5,1,5)
plot(T(1:(end-1)), r_out)
hold on
plot([0 T(end)], [1 1], 'r')
ylim([-0.2 1.2])
xlabel('t')
ylabel('r')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'7_3_TestEquation','-dpdf','-r0')

%% 1.4.Evolution of the error vs. stepsize
clc
close all
tspan = [0 30];
tols = logspace(-9, -4, 15);
h0 = 0.1;

e_ls_ex = zeros(1,size(tols,2));
e_gs_ex = zeros(1,size(tols,2));
e_ls_im = zeros(1,size(tols,2));
e_gs_im = zeros(1,size(tols,2));

i = 0;
for tol=tols
    i = i+1;
    [T,X,X_real,e_l,e_g,e_out,r_out,h_out,info] = DOPRI54TE_adaptive(@TestEquation, ...
    @Exponential, tspan, h0, x0,tol,tol, butcher, args);
    
    e_ls_ex(i) = mean(e_l);
    e_gs_ex(i) = mean(e_g);
end

figure()
loglog(tols,e_ls_ex)
hold on
loglog(tols, 5*tols.^(1.15),'LineStyle','--','Color','k')
xlabel('tol')
ylabel('local error')
legend('DOPRI54','O(tol^{1.15})','Location','NorthWest')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'7_3_localerror','-dpdf','-r0')

figure()
loglog(tols,e_gs_ex)
hold on
loglog(tols, 1*tols.^(1),'black','LineStyle','--','Color','k')
xlabel('tol')
ylabel('global error')
legend('DOPRI54','O(tol)','Location','NorthWest')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'7_3_globalerror','-dpdf','-r0')

%% Stability
clc
clear 
close all
method = 'DOPRI54';
alpha = -5:0.01:5;
beta =  -5:0.01:5;

solver = ERKSolverErrorEstimationParameters(method);
A = solver.AT';
b = solver.b;
c = solver.d;
d = solver.d;

nreal = length(alpha);
nimag = length(beta);
I = eye(size(A));
e = ones(size(A,1),1);

for kreal = 1:nreal
    for kimag = 1:nimag
        z = alpha(kreal) + i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b'*tmp;
%         Ehat = z*d'*tmp;
        f = exp(z);
        E = R-f;
%         EhatmE = Ehat-E;
        absR(kimag,kreal) = abs(R);
%         absEhatmE(kimag,kreal) = abs(EhatmE);
%         absEhat(kimag,kreal)   = abs(Ehat);
        absE(kimag,kreal) = abs(E);
        absF(kimag,kreal) = abs(f);
    end
end

figure()
imagesc(alpha,beta,absR,[0 1]);
hold on
grid on
colorbar
colormap jet
axis image
axis xy
plot([alpha(1) alpha(end)], [0 0], 'LineWidth',1,'Color','black')
plot([0 0], [alpha(1) alpha(end)], 'LineWidth',1,'Color','black')
xticks(-5:1:5);
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'7_3_stability_regions','-dpdf','-r0')

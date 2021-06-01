%% 1.3.Test Equation error measurements
clc
clear
close all

tspan = [0 5];
h = 0.2;
x0 = 1;
lambda = -1;

args = {lambda};

[T,X,X_real,e_l,e_g] = ClassicRK4TE_fixed(@TestEquation, @Exponential, tspan, h, x0, args);

figure()
yyaxis left
plot(T,X)
hold on
plot(T,X_real,'--','Color',[0.9290, 0.6940, 0.1250])

xlabel('t')
ylabel('x(t)')

yyaxis right
plot(T,e_l)
plot(T,e_g)
hold off
ylabel('errors')
legend('numerical solution','real solution', 'local error', 'global error')
limits = ylim;

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'5_3_TestEquation','-dpdf','-r0')

%% 1.4.Evolution of the error vs. stepsize
clc
close all
tspan = [0 30];
hs = 0.01:0.005:0.5;

e_ls_ex = zeros(1,size(hs,2));
e_gs_ex = zeros(1,size(hs,2));
e_ls_im = zeros(1,size(hs,2));
e_gs_im = zeros(1,size(hs,2));

i = 0;
for h=hs
    i = i+1;
    [T,X,X_real,e_l,e_g] = ClassicRK4TE_fixed(@TestEquation, @Exponential, tspan, h, x0, args);
    
    e_ls_ex(i) = mean(e_l);
    e_gs_ex(i) = mean(e_g);
end

figure()
loglog(hs,e_ls_ex)
hold on
loglog(hs, 0.025*hs.^5,'LineStyle','--','Color','k')
xlabel('h')
ylabel('local error')
legend('RK4','O(h^5)','Location','NorthWest')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'5_3_localerror','-dpdf','-r0')

figure()
loglog(hs,e_gs_ex)
hold on
loglog(hs, 0.03*hs.^4,'black','LineStyle','--','Color','k')
xlabel('h')
ylabel('global error')
legend('RK4','O(h^4)','Location','NorthWest')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'5_3_globalerror','-dpdf','-r0')

%% Stability
clc
clear 
close all
method = 'RK4C';
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
print(gcf,'5_3_stability_regions','-dpdf','-r0')

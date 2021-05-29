%% 1.3.Test Equation error measurements
clc
close all

tspan = [0 5];
h = 0.2;
x0 = 1;
lambda = -1;

args = {lambda};

[T,X,X_real,e_l,e_g] = EulerExplicit(@TestEquation, @Exponential, tspan, h, x0, args);
[T1,X1,X_real1,e_l1,e_g1] = EulerImplicitTest(@Exponential, tspan, h, x0, lambda);

figure()
yyaxis left
plot(T,X)
hold on
plot(T,X_real)

xlabel('t')
ylabel('x(t)')

yyaxis right
plot(T,e_l)
plot(T,e_g)
hold off
ylabel('errors')
legend('numerical solution','real solution', 'local error', 'global error')
limits = ylim;

% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters 1_explicit_0_2

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'1_explicit_0_2','-dpdf','-r0')

figure()
yyaxis left
plot(T1,X1)
hold on
plot(T1,X_real1)

xlabel('t')
ylabel('x(t)')

yyaxis right
plot(T1,e_l1)
plot(T1,e_g1)
ylim(limits)
hold off
ylabel('errors')
legend('numerical solution','real solution', 'local error', 'global error')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'1_implicit_0_2','-dpdf','-r0')


%% 1.4.Evolution of the error vs. stepsize
clc
close all
tspan = [0 20];
hs = 0.001:0.005:0.5;

e_ls_ex = zeros(1,size(hs,2));
e_gs_ex = zeros(1,size(hs,2));
e_ls_im = zeros(1,size(hs,2));
e_gs_im = zeros(1,size(hs,2));

i = 0;
for h=hs
    i = i+1;
    [T,X,X_real,e_l,e_g] = EulerExplicit(@TestEquation, @Exponential, tspan, h, x0, args);
    [T1,X1,X_real1,e_l1,e_g1] = EulerImplicitTest(@Exponential, tspan, h, x0, lambda);
    
    idx = find(abs(T-t_plot) == min(abs(T-t_plot)));
    idx = idx(1);
    e_ls_ex(i) = mean(e_l);
    e_gs_ex(i) = mean(e_g);
    
    e_ls_im(i) = mean(e_l1);
    e_gs_im(i) = mean(e_g1);
end

figure()
plot(hs,e_ls_ex)
hold on
plot(hs,e_ls_im)
plot(hs, 0.025*hs.^2,'LineStyle','--','Color','k')
% ax = gca;
% ax.ColorOrderIndex = 1;
xlabel('h')
ylabel('local error')
legend('explicit', 'implicit','O(h^2)','Location','NorthWest')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'1_4_localerror','-dpdf','-r0')

figure()
plot(hs,e_gs_ex)
hold on
% ax.ColorOrderIndex = 2;
plot(hs,e_gs_im)
plot(hs, 0.03*hs,'black','LineStyle','--','Color','k')
xlabel('h')
ylabel('global error')
legend('explicit', 'implicit','O(h)','Location','NorthWest')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'1_5_globalerror','-dpdf','-r0')

%% 1.5.Stability Explicit and Implicit
clc
clear
close all

real = linspace(-4,1,1e3);
img = linspace(-2,2,1e3);
[X,Y] = meshgrid(real,img);

Z = -sqrt((X+1).^2 + (Y).^2); 

figure()
contourf(X,Y,Z,[-1 -1]);
colormap winter
hold on
plot(-1, 0, ".k");
plot(axis, axis.*0, "k");
plot(axis.*0, axis, "k");
xlim([-4 1]);
ylim([-2 2]);
xlabel('Re(h\lambda)')
ylabel('Im(h\lambda)')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'1_6_explicit','-dpdf','-r0')

% Implicit
f = @(z) 1./(1-z);
real = linspace(-1,4,1e3);
img = linspace(-2,2,1e3);
[X,Y] = meshgrid(real,img);

Z = -abs(f(X + Y*i)); 

figure()
contourf(X,Y,Z,[-1 -1]);
colormap winter
hold on
plot(1, 0, ".k");
plot(axis, axis.*0, "k");
plot(axis.*0, axis, "k");
xlim([-1 4]);
ylim([-2 2]);
xlabel('Re(h\lambda)')
ylabel('Im(h\lambda)')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'1_6_implicit','-dpdf','-r0')
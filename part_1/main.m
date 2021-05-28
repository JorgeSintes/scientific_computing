%% Test Equation error measurements
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


%% Evolution of the error vs. stepsize
clc
close all
tspan = [0 15];
hs = 0.005:0.001:0.5;

e_ls_ex = zeros(1,size(hs,2));
e_gs_ex = zeros(1,size(hs,2));
e_ls_im = zeros(1,size(hs,2));
e_gs_im = zeros(1,size(hs,2));

i = 0;
for h=hs
    i = i+1;
    [T,X,X_real,e_l,e_g] = EulerExplicit(@TestEquation, @Exponential, tspan, h, x0, args);
    [T1,X1,X_real1,e_l1,e_g1] = EulerImplicitTest(@Exponential, tspan, h, x0, lambda);
    e_ls_ex(i) = sum(e_l);
    e_gs_ex(i) = sum(e_g)/size(e_g,2);
    e_ls_im(i) = sum(e_l1);
    e_gs_im(i) = sum(e_g1)/size(e_g1,2);
end

figure()
plot(hs,e_ls_ex)
hold on
plot(hs,e_ls_im)
% ax = gca;
% ax.ColorOrderIndex = 1;
legend('explicit', 'implicit')

figure()
plot(hs,e_gs_ex, '-.')
hold on
% ax.ColorOrderIndex = 2;
plot(hs,e_gs_im, '-.')

xlabel('h')
ylabel('error')
legend('explicit', 'implicit')

%% Trial
clc
close all
h = hs(1);
h = 0.00001;

[T,X,X_real,e_l,e_g] = EulerExplicit(@TestEquation, @Exponential, tspan, h, x0, args);
[T1,X1,X_real1,e_l1,e_g1] = EulerImplicitTest(@Exponential, tspan, h, x0, lambda);

figure()
plot(T,X)
hold on
plot(T,X1)
plot(T,X_real)
legend('explicit', 'implicit', 'real')
sum(e_l)
sum(e_g)
sum(e_l1)
sum(e_g1)
sum(X-X_real)
sum(X1-X_real1)

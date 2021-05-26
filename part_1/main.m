%% Test Equation error measurements
clc
close all

tspan = [0 5];
h = 0.1;
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
legend('real solution', 'numerical solution', 'local error', 'global error')
limits = ylim;

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
legend('real solution', 'numerical solution', 'local error', 'global error')


%% Evolution of the error vs. stepsize
tspan = [0 15];
hs = 0.005:0.005:1;

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
    e_gs_ex(i) = sum(e_g);
    e_ls_im(i) = sum(e_l1);
    e_gs_im(i) = sum(e_g1);
end

figure()
plot(hs,e_ls_ex)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(hs,e_gs_ex, '-.')

plot(hs,e_ls_im)
ax.ColorOrderIndex = 2;
plot(hs,e_gs_im, '-.')

xlabel('h')
ylabel('error')
legend('explicit - local', 'explicit - global', 'implicit - local', 'implicit - global')

%% Comparison Van der Pol mu=15
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
tols =  [1e-2 1e-3 1e-4 1e-5 1e-6];

figure('Renderer', 'painters', 'Position', [10 10 1000 800])

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    [T,X,r_out,h_out,info] = EulerExplicit_adaptive(@VanderPol, tspan, h0, x0, abstol, reltol, args);
    infos(:,count) = info;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on
xlabel('tol');
ylabel('Function Evaluations');
set(gca,'Xdir','reverse');

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on
xlabel('tol');
ylabel('Calculated steps');
set(gca,'Xdir','reverse');

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on
xlabel('tol');
ylabel('Accepted steps');
set(gca,'Xdir','reverse');

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on
xlabel('tol');
ylabel('Rejected steps');
set(gca,'Xdir','reverse');


%% Implicit
infos = zeros(4,size(tols,2));
count = 0;

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    [T,X,r_out,h_out,info] = EulerImplicit_adaptive(@VanderPol_fJ, tspan, h0, x0, abstol, reltol, args);
    infos(:,count) = info;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% Classical Runge-Kutta

count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    [T,X,r_out,h_out,info] = ClassicRK4_adaptive(@VanderPol, tspan, h0, x0, abstol, reltol, args);
    infos(:,count) = info;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% DOPRI54

butcher = ERKSolverErrorEstimationParameters('DOPRI54');
count = 0;
infos = zeros(4,size(tols,2));

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    [T,X,e_out,r_out,h_out,info] = DOPRI54_adaptive(@VanderPol, tspan, h0, x0, abstol, reltol, butcher, args);
    infos(:,count) = info;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% ESDIRK23

method = 'ESDIRK23';
fun  = 'VanderPolFun';
jac  = 'VanderPolJac';
t0 = tspan(1);
tf = tspan(end);

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
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% ode45
count = 0;
infos = zeros(4,size(tols,2));
info = zeros(4,1);

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
    sol = ode45(@(t,x) VanderPol(t,x,mu),tspan, x0, opts);
    info(1) = sol.stats.nfevals;
    info(3) = sol.stats.nsteps;
    info(4) = sol.stats.nfailed;
    info(2) = info(3) + info(4);
    infos(:,count) = info;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on

%% ode15s
count = 0;
infos = zeros(4,size(tols,2));
info = zeros(4,1);

for tol=tols
    count = count + 1;
    abstol = tol;
    reltol = tol;

    opts = odeset('RelTol',tol,'AbsTol',tol,'InitialStep',h0);
    sol = ode15s(@(t,x) VanderPol(t,x,mu),tspan, x0, opts);
    info(1) = sol.stats.nfevals;
    info(3) = sol.stats.nsteps;
    info(4) = sol.stats.nfailed;
    info(2) = info(3) + info(4);
    infos(:,count) = info;
end

subplot(2,2,1)
loglog(tols,infos(1,:),'o-');
hold on

subplot(2,2,2)
loglog(tols,infos(2,:),'o-');
hold on

subplot(2,2,3)
loglog(tols,infos(3,:),'o-');
hold on

subplot(2,2,4)
loglog(tols,infos(4,:),'o-');
hold on
legend('Ex. Euler', 'Im. Euler','RK4','DOPRI54','ESDIRK23','ode45','ode15s','Location','SouthWest');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'9_mu_15','-dpdf','-r0')
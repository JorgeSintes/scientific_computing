%% Van der Pol Euler-Maruyama State independent
clc
clear
close all

Ns = 5;
sigma_EM_in = 1;
sigma_IE_in = 1;
sigma_EM_dep = 0.5;
sigma_IE_dep = 0.5;

sigma = sigma_EM_in;


tspan = [0 40];
h = 0.01;
x0 = [0.5;0.5];
nW = 1;
mu = 3;

args = {mu, sigma};

for i=1:Ns
    [T,X(:,:,i)] = EulerMaruyama(@VanderPolDrift,@VanderPolDiffusion1,tspan,h,x0,nW,args);
end

[T, Xd] = EulerMaruyama(@VanderPolDrift,@VanderPolDiffusion1,tspan,h,x0,nW,{mu,0});

figure('Renderer', 'painters', 'Position', [10 10 1200 533])

for i=1:Ns
    subplot(2,2,1)
    plot(T,X(1,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,3)
    plot(T,X(2,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,[2,4])
    plot(X(1,:,i),X(2,:,i))
    hold on
    
    xlabel('x_1');
    ylabel('x_2');
end

subplot(2,2,1)
plot(T,Xd(1,:),'LineWidth',2,'Color','black');

subplot(2,2,3)
plot(T,Xd(2,:),'LineWidth',2,'Color','black');

subplot(2,2,[2,4])
plot(Xd(1,:),Xd(2,:),'LineWidth',2,'Color','black')
XL = get(gca, 'XLim');
maxlim = max(abs(XL));
set(gca, 'XLim', [-maxlim maxlim]);
YL = get(gca, 'YLim');
maxlim = max(abs(YL));
set(gca, 'YLim', [-maxlim maxlim]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_5_EulerMaruyama_independent','-dpdf','-r0')

%% Van der Pol Implicit-Explicit State Independent
clc
close all

sigma = sigma_IE_in;
args = {mu, sigma};

for i=1:Ns
    [T,X(:,:,i)] = ImplicitExplicit(@VanderPolDrift,@VanderPolDiffusion1,tspan,h,x0,nW,args);
end

[T, Xd] = ImplicitExplicit(@VanderPolDrift,@VanderPolDiffusion1,tspan,h,x0,nW,{mu,0});

figure('Renderer', 'painters', 'Position', [10 10 1200 533])

for i=1:Ns
    subplot(2,2,1)
    plot(T,X(1,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,3)
    plot(T,X(2,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,[2,4])
    plot(X(1,:,i),X(2,:,i))
    hold on
    
    xlabel('x_1');
    ylabel('x_2');
end

subplot(2,2,1)
plot(T,Xd(1,:),'LineWidth',2,'Color','black');

subplot(2,2,3)
plot(T,Xd(2,:),'LineWidth',2,'Color','black');

subplot(2,2,[2,4])
plot(Xd(1,:),Xd(2,:),'LineWidth',2,'Color','black')
XL = get(gca, 'XLim');
maxlim = max(abs(XL));
set(gca, 'XLim', [-maxlim maxlim]);
YL = get(gca, 'YLim');
maxlim = max(abs(YL));
set(gca, 'YLim', [-maxlim maxlim]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_5_ImplicitExplicit_independent','-dpdf','-r0')

%% Van der Pol Euler-Maruyama State Dependent
clc
close all

sigma = sigma_EM_dep;
args = {mu, sigma};

for i=1:Ns
    [T,X(:,:,i)] = EulerMaruyama(@VanderPolDrift,@VanderPolDiffusion2,tspan,h,x0,nW,args);
end

[T, Xd] = EulerMaruyama(@VanderPolDrift,@VanderPolDiffusion2,tspan,h,x0,nW,{mu,0});

figure('Renderer', 'painters', 'Position', [10 10 1200 533])

for i=1:Ns
    subplot(2,2,1)
    plot(T,X(1,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,3)
    plot(T,X(2,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,[2,4])
    plot(X(1,:,i),X(2,:,i))
    hold on
    
    xlabel('x_1');
    ylabel('x_2');
end

subplot(2,2,1)
plot(T,Xd(1,:),'LineWidth',2,'Color','black');

subplot(2,2,3)
plot(T,Xd(2,:),'LineWidth',2,'Color','black');

subplot(2,2,[2,4])
plot(Xd(1,:),Xd(2,:),'LineWidth',2,'Color','black')
XL = get(gca, 'XLim');
maxlim = max(abs(XL));
set(gca, 'XLim', [-maxlim maxlim]);
YL = get(gca, 'YLim');
maxlim = max(abs(YL));
set(gca, 'YLim', [-maxlim maxlim]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_5_EulerMaruyama_dependent','-dpdf','-r0')

%% Van der Pol Implicit-Explicit State dependent
clc
close all

sigma = sigma_IE_dep;
args = {mu, sigma};

for i=1:Ns
    [T,X(:,:,i)] = ImplicitExplicit(@VanderPolDrift,@VanderPolDiffusion2,tspan,h,x0,nW,args);
end

[T, Xd] = ImplicitExplicit(@VanderPolDrift,@VanderPolDiffusion2,tspan,h,x0,nW,{mu,0});

figure('Renderer', 'painters', 'Position', [10 10 1200 533])

for i=1:Ns
    subplot(2,2,1)
    plot(T,X(1,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,3)
    plot(T,X(2,:,i));
    hold on
    xlabel('t');
    ylabel('x_2');
    
    subplot(2,2,[2,4])
    plot(X(1,:,i),X(2,:,i))
    hold on
    
    xlabel('x_1');
    ylabel('x_2');
end

subplot(2,2,1)
plot(T,Xd(1,:),'LineWidth',2,'Color','black');

subplot(2,2,3)
plot(T,Xd(2,:),'LineWidth',2,'Color','black');

subplot(2,2,[2,4])
plot(Xd(1,:),Xd(2,:),'LineWidth',2,'Color','black')
XL = get(gca, 'XLim');
maxlim = max(abs(XL));
set(gca, 'XLim', [-maxlim maxlim]);
YL = get(gca, 'YLim');
maxlim = max(abs(YL));
set(gca, 'YLim', [-maxlim maxlim]);

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'4_5_ImplicitExplicit_dependent','-dpdf','-r0')

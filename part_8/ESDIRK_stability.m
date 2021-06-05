%% Stability
clc
clear 
close all
method = 'ESDIRK23';
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
print(gcf,'8_2_stability_regions','-dpdf','-r0')

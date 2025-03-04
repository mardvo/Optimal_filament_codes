close all 
clear all


%%% 
load('pw_sqrt_span_a1_V1.mat')

Ff = NaN*zeros([length(a1s)-3,length(Volumes)]);
X0 = NaN*zeros([length(a1s)-3,length(Volumes)]);

for s=1:length(a1s)-3
    for j=1:length(Volumes)
    
        [Ff(s,j),ind1]=(max(-F_prop(s,j,:)));
        X0(s,j) = x(num(ind1));
    
    end
end

figure
[X1,X2]=meshgrid(a1s(1:end-3).^4,Volumes);
contourf(X2,X1,Ff',100,'Edgecolor','none')
cc=colorbar;
ylabel('$a$','Interpreter','Latex')
xlabel('$V_0$','Interpreter','Latex')
set(gca,'Fontsize',16)
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel(cc, '$F$','Interpreter','Latex')
xlim([min(Volumes),max(Volumes)])
ylim([10^(-3),max(a1s(1:end-3).^4)])

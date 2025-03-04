close all 
clear all

load('pw_const_a1a2.mat') 
Ff = NaN*zeros([length(a1s),length(a2s)]);
X0 = NaN*zeros([length(a1s),length(a2s)]);
for m=1:length(a2s)
    for mm=1:m-1
    [Ff(mm,m),ind1]=(max(abs(F(mm,m,:))));
    X0(mm,m) = x0(ind1);
    dF(mm,m,:)=first_der(F(mm,m,:),x0);
    end
end

% lower boundaries is an array a1s, 
% upper boundaries is an array a2s

% figure 5a
 [X1,X2]=meshgrid(a1s,a2s);
contourf(X2,X1,Ff',100,'EdgeColor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
c=colorbar;
ylabel('$a$','Interpreter','Latex')
xlabel('$b$','Interpreter','Latex')
set(gca,'Fontsize', 16)
xlim([10^(-3),10])
ylim([10^(-3),10])
ylabel(c, '$F$','Interpreter','Latex','FontSize',18)

%print(gcf,'F_max.png','-dpng','-r600'); 

% figure 5b
figure 
contourf(X2,X1,X0',100,'EdgeColor','none')
set(gca,'XScale','log')
set(gca,'YScale','log')
c=colorbar;
ylabel(c, '$x_0$','Interpreter','Latex','FontSize',18)
ylabel('$a$','Interpreter','Latex')
xlabel('$b$','Interpreter','Latex')
xlim([10^(-3),10])
ylim([10^(-3),10])
set(gca,'Fontsize', 16)
%print(gcf,'x0.png','-dpng','-r600'); 



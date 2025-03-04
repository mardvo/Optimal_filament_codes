clear all;
close all;

npts=1000;

ni=150;
a1=1e-3;
a2=1; 

x = linspace(0,1,500);

nx0=200;
x0=linspace(5/nx0,1-5/nx0,nx0);
F=zeros([ni,nx0]);

for i=1:length(x0)
    [~,~,F(i),~,~,F1(i),F2(i)] = const_2_analit(x0(i),npts,a1,a2);

    [~,~,F_end(i),~,~] = const1_analit(x0(i),x,npts,a1/(1-x0(i))^4);
end


figure
plot(x0,-F1,'Linewidth',2)
hold on 
plot(x0,-F2,'Linewidth',2)
plot(x0,-F_end.*(1-x0).^2,'k--','Linewidth',4)
xline(0.6742)

ll=legend('$F_1$','$F_2$','$F_{end}$');
ll.Interpreter = 'Latex';
ll.FontSize = 20; 
xlabel('$x_0$','Interpreter','Latex')
%ylabel('$F$','Interpreter','Latex')
set(gca,'Fontsize', 20)


clear all;
close all;

npts=1000;

x = linspace(0,1,npts);
k=0;
A0s=10.^linspace(-4,2,200);

for A0=A0s
k=k+1;
consts; % solve main problem
const_cong; % solve the conjugate problem
Fp(k) = F_prop;
end

figure
semilogx(A0s,-Fp,'b','Linewidth',3)
hold on

ylabel('$F$','Interpreter','Latex')

xlabel('$A$','Interpreter','Latex')
%ylabel('$dF/|dF|$','Interpreter','Latex')
set(gca,'Fontsize',16)

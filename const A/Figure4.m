clear all;
close all;

npts=1000;

x = linspace(0,1,npts);
k=0;

linetypes = ["b-","b-.","b--","b:"];

for A0=[1e-3, 1e-2, 1e-1]
k=k+1;
consts; % solve main problem
const_cong; % solve the conjugate problem
dF1(k,:)=dFdA; 
dF_norm(k)= sqrt(trapz(x,dFdA.^2))
end

figure
for k=1:3
plot(x,-dF1(k,:)'./dF_norm(k)',linetypes(k),'Linewidth',2)
hold on

end
ylabel('$(dF/dA)/|dF/dA|$','Interpreter','Latex')
ylim([-1.5,3])

plot([0,1],[0,0],'k-')
ll = legend('$A=10^{-3}$','$A=10^{-2}$','$A=10^{-1}$');

ll.Interpreter = 'Latex';

xlabel('$x$','Interpreter','Latex')
%ylabel('$dF/|dF|$','Interpreter','Latex')
set(gca,'Fontsize',16)

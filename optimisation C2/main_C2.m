%%% elastic flagellum - first attempt to code it --- bvp4c
% everything is dimensionless.
% with Fourier modes

clear all
close all

global l a1

l=1; 
a1=(0.01)^(1/4); 
V0=0.5; % fixed volume
rho0=1e-2; % initial penalty parameter

npts=500;
x = linspace(0,l,npts);

rin=A_cubic(a1,1,x); % initial condition

[Rr,F_prop,dF]=descent_penalty(rin,V0,rho0);

% 
figure
%plot(x,Rr.^4')
hold on 
plot(x,Rr(1,:).^4','k--')
plot(x,Rr(end,:).^4','k','LineWidth',3)
xlabel('x')
ylabel('A')
set(gca,'Fontsize',14)

save('opt_C2.mat')
function [r_rho,F_rho,dF] = descent_penalty(rin,V,rho)
global l a1
npts=500;
x = linspace(0,l,npts);

incon = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
solinit = bvpinit(linspace(0,l,npts),incon);
options = bvpset('RelTol',1e-5,'Stats','on','NMax',10000);
gamma=0.001;

    r(1,:)=rin/sqrt(trapz(x,rin.^2))*sqrt(V);%rescale the first point to match the desired volume
    
    xx=x;

    sol = bvp5c(@(z,y)twoode(z,y,interp1(xx,r(1,:).^4,z)),@twobc,solinit,options); % solve ODEs
    y2 = deval(sol,x);
    F_prop(1,:) = trapz(x,y2(2,:).*y2(5,:)-y2(6,:).*y2(1,:))/2+rho*trapz(x,r(1,:).^2-V)^2; % compute the objective function
    dF(1,:)=4.*(y2(15,:).*y2(3,:)+y2(11,:).*y2(7,:))./r(1,:).^5+4*rho*trapz(x,r(1,:).^2-V)*r(1,:); % dF/dr
 
s=1;
err3=1;
delta3=5e-5;
delta2 = 1e-2;
delta1 = 1e-2;
F_rho(1)=F_prop(1,:); r_rho(1,:)=r(1,:); 
while ((err3>delta3 || s<12) & s<20)
    rin=r_rho(s,:);
    s=s+1

    delta2=max(delta2/2,1e-4);% reduce the error
    delta1=max(delta1/4,1e-4);% reduce the error
    
    [r,F_prop,dF]=proj_armijo(rin,V,rho,delta1,delta2);

F_rho(s)=F_prop(end,:);
r_rho(s,:)=r(end,:);
err3=abs((F_rho(s)-F_rho(s-1))/F_rho(s-1))+max(abs(r_rho(s,:)-r_rho(s-1,:))) % calculate err3
rho=rho*2; % increase the penalty parameter
end
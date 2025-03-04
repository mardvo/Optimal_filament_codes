%%% elastic flagellum - first attempt to code it --- bvp4c
% everything is dimensionless.
clear all
close all

global l a1 a2

l=1; 
a2=0.1; % upper boundary (=b in the paper)
a1=0.01; % lower boundary (=a in the paper)

npts=500; % points in x
toll=1e-6; % tolerance for ODEs
x = linspace(0,l,npts);

incon = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
solinit = bvpinit(linspace(0,l,npts),incon);
options = bvpset('RelTol',1e-8,'Stats','on','NMax',10000);
y2_p = solinit.y+1;
y2 = solinit.y;
n=0;

gamma=0.001;

delta = 1e-5;
k=1;
ni=500;
    A1=A_cubic(a1,a2,x); 
    A(1,:)=A1;% initial condition
    xx=x;
    solinit.y = y2_p;
    sol = bvp5c(@(z,y)twoode(z,y,interp1(xx,A(1,:),z)),@twobc,solinit,options); % solve odes
    y2 = deval(sol,x);
    F_prop(1,:) = trapz(x,y2(2,:).*y2(5,:)-y2(6,:).*y2(1,:))/2 % calculate propulsive force
    dF(1,:)=(y2(15,:).*y2(3,:)+y2(11,:).*y2(7,:))./A(1,:).^2; % calculate dF/dA   
    norm_dF(1)=sqrt(trapz(x,dF(1,:).^2));% calculate |dF/dA|
k=1;
err1=1;
err2=1;
while (err1>delta & err2>delta)
    % choose the time step
    k=k+1;
    sigma1=1;
    sigma=sigma1;
    Fp=F_prop(k-1,:);
    % wk
    A(k,:)=Proj(A(k-1,:)-sigma*dF(k-1,:)/norm_dF(k-1),a1,a2);
    
    while  Fp-F_prop(k-1,:)>-gamma/sigma*trapz(x,(A(k,:)-A(k-1,:)).^2) 
    sigma=sigma1;
    A(k,:)=Proj(A(k-1,:)-sigma*dF(k-1,:)/norm_dF(k-1),a1,a2);     

        try
            k
    solinit.y = y2;
    sol = bvp5c(@(z,y)twoode(z,y,interp1(xx,A(k,:),z)),@twobc,solinit,options);
    y2 = deval(sol,x);
    Fp = trapz(x,y2(2,:).*y2(5,:)-y2(6,:).*y2(1,:))/2
        catch
        end
    sigma1=sigma1/2;
    end
    %A(k,:)=A(k-1,:)-sigma*dF(k-1,:)./norm_dA(k-1); 

    F_prop(k,:) = Fp;
    dF(k,:)=(y2(15,:).*y2(3,:)+y2(11,:).*y2(7,:))./A(k,:).^2;
    norm_dF(k)=sqrt(trapz(x,dF(k,:).^2));    
    err2=abs((F_prop(k,:)-F_prop(k-1,:))/F_prop(k,:));

    A(k,:)=Proj(A(k-1,:)-sigma*dF(k-1,:)/norm_dF(k-1),a1,a2);
    err1=(max(abs(A(k,:)-A(k-1,:))));
end

figure
plot(x,A')
xlabel('x')
ylabel('A')
set(gca,'Fontsize',14)

save opt_reg1.mat

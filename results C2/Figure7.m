close all 
clear all

%
load("rand_001_0.2.mat")
%nv=5;
Volumes = a1^2*[2:5,7,9,10];
for nv=1:6
    kkk=1;
    filename="rand_001_"+num2str(Volumes(nv))+".mat";
    load(filename)
    %Rrin(kkk,:)=(Rr(1,:),a1,a2);
    Rend(nv,:)=Rr(end,:);
    F(nv)=F_prop(end)
    %Volumes = (a2^2-a1^2)*[1:9]/10+a1^2;
end

figure
plot(x,Rend.^4','Linewidth',2)
hold on 
%plot(x,Rin','--')
xlabel('$x$','Interpreter','Latex')
ylabel('$A$','Interpreter','Latex')
set(gca,'Fontsize',16)

figure
for k=1:nv
plot(Volumes(k),-F(k),'o','Linewidth',6)
hold on
end
xlabel('$V_0$','Interpreter','Latex')
ylabel('$F$','Interpreter','Latex')

set(gca,'Fontsize',16)


load("rand_001_0.5.mat")


incon = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
solinit = bvpinit(linspace(0,1,npts),incon);
options = bvpset('RelTol',1e-8,'Stats','on','NMax',10000);
    sol = bvp5c(@(z,y)twoode(z,y,interp1(x,Rr(end,:).^4,z)),@twobc,solinit,options);
    y2 = deval(sol,x);
    omega=1;
    t=linspace(0,2*pi,15);
%t=[0,pi/2,pi,3*pi/2]/omega
% Create a grayscale colormap
colormap(gray(256));

% Generate a vector of grayscale values for line colors
gray_colors = linspace(0.1, 0.9, 14);

figure
for i=2:14
    Y(:,i)=y2(1,:)*cos(omega*t(i))+y2(5,:)*sin(omega*t(i));
    line_color = [gray_colors(i) gray_colors(i) gray_colors(i)];
    plot(x, Y(:, i), 'Color', line_color, 'LineWidth', 2);
    hold on;

end
Y(:,1)=y2(1,:)*cos(omega*t(1))+y2(5,:)*sin(omega*t(1));
set(gca,'Fontsize',30)
xlabel('$x$','Interpreter','Latex','Fontsize',30)
ylabel('$y$','Interpreter','Latex','Fontsize',30)
plot(x, Y(:, 1), 'k-','LineWidth', 3);
% Adjust the colormap range to grayscale
colormap gray;

% Hold off to prevent further plots from overwriting the current one
hold off;


%%% 
load('pw_sqrt_001.mat')

[Ff,ind1]=(max(-F_prop'));
[Fff,ind2]=max(Ff)
Volumes = a1^2*[2:5,7,9,10];
for nv=1:6
A1(nv,:)=r_pw_sqrt(a1,Volumes(nv),x,x(num(ind2+1))); 
incon = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
solinit = bvpinit(linspace(0,1,npts),incon);
options = bvpset('RelTol',1e-8,'Stats','on','NMax',10000);
    sol = bvp5c(@(z,y)twoode(z,y,interp1(x,A1(nv,:).^4,z)),@twobc,solinit,options);
    y2 = deval(sol,x);
    F_quad(nv) = trapz(x,y2(2,:).*y2(5,:)-y2(6,:).*y2(1,:))/2;
end

figure
plot(x,A1,'--')
hold on
plot(x,Rend','Linewidth',2)
hold on 
%plot(x,Rin','--')
xlabel('x')
ylabel('r')
set(gca,'Fontsize',14)


figure
plot(x,A1(4,:),'k--','Linewidth',6)
hold on
plot(x,Rend(4,:)','Linewidth',3)
hold on 

xlabel('$x$','Interpreter','Latex')
ylabel('$r$','Interpreter','Latex')
set(gca,'Fontsize',32)
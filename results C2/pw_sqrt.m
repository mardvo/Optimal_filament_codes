%%% span an approximated shape for a1, V0 and x0

clear all
close all

l=1;


npts=500;
toll=1e-6;
x = linspace(0,l,npts);

incon = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
solinit = bvpinit(linspace(0,l,npts),incon);
options = bvpset('RelTol',1e-8,'Stats','on','NMax',10000);
y2_p = solinit.y+1;
y2 = solinit.y;
n=0;
k=1;
ni=50;

num=linspace(10,npts,ni);
a1s=(10.^[-3:0.02:-1]).^(1/4);
Volumes=0.001^(1/2)*10.^[0.1:0.02:1.7];
F_prop=zeros([length(a1s),length(Volumes),ni]);
for s=1:length(a1s)
    a1=a1s(s); 
    for j=1:length(Volumes)
            V=Volumes(j);
    
        for i=1:length(num)
            if V<=a1^2
                F_prop(s,j,i)=NaN;
            else
                A1=r_pw_sqrt(a1,V,x,x(num(i)));    
                A(i,:)=A1;
                xx=x;
                solinit.y = y2_p;
                sol = bvp5c(@(z,y)twoode(z,y,interp1(xx,A(i,:).^4,z)),@twobc,solinit,options);
                y2 = deval(sol,x);
                F_prop(s,j,i) = trapz(x,y2(2,:).*y2(5,:)-y2(6,:).*y2(1,:))/2;
                dF(i,:)=(y2(15,:).*y2(3,:)+y2(11,:).*y2(7,:))./A(i,:).^2;
                norm_dA(i,j,s)=sqrt(trapz(x,dF(i,:).^2));
            end
        end
    end
end

save pw_sqrt_span_a1_V1.mat


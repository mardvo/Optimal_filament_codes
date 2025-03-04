function [r,F_prop,dF] = proj_armijo(rin,V,rho,delta1,delta2)

global l a1

npts=500;
x = linspace(0,l,npts);

incon = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
solinit = bvpinit(linspace(0,l,npts),incon);
options = bvpset('RelTol',1e-8,'Stats','on','NMax',1000);
gamma=0.001;

    r(1,:)=rin;
    xx=x;

    sol = bvp5c(@(z,y)twoode(z,y,interp1(xx,r(1,:).^4,z)),@twobc,solinit,options);
    y2 = deval(sol,x);
    F_prop(1,:) = trapz(x,y2(2,:).*y2(5,:)-y2(6,:).*y2(1,:))/2+rho*trapz(x,r(1,:).^2-V)^2;
    dF(1,:)=4.*(y2(15,:).*y2(3,:)+y2(11,:).*y2(7,:))./r(1,:).^5+4*rho*trapz(x,r(1,:).^2-V)*r(1,:); % dF/dr

    norm_dF(1)=sqrt(trapz(x,dF(1,:).^2));

k=1;
err1=1;
err2=1; err3=1;
while ((err1>delta1 || err2>delta2) & k<1000)
    % choose the step size
    k=k+1;
    sigma1=1;
    sigma=sigma1;
    Fp=F_prop(k-1,:);
    % wk
    r(k,:)=Proj(r(k-1,:)-sigma*dF(k-1,:)/norm_dF(k-1),a1);
    
    while  Fp-F_prop(k-1,:)>-gamma/sigma*trapz(x,(r(k,:)-r(k-1,:)).^2) 
    sigma=sigma1;
    r(k,:)=Proj(r(k-1,:)-sigma*dF(k-1,:)/norm_dF(k-1),a1);     

        try
            k
    
    solinit.y = y2;
    sol = bvp5c(@(z,y)twoode(z,y,interp1(xx,r(k,:).^4,z)),@twobc,solinit,options);
    y2 = deval(sol,x);
    Fp = trapz(x,y2(2,:).*y2(5,:)-y2(6,:).*y2(1,:))/2+rho*trapz(x,r(k,:).^2-V)^2
        catch
        end
    sigma1=sigma1/2;
    if sigma1<1e-6
        return
    end
    end
    F_prop(k,:) = Fp;
    dF(k,:)=4.*(y2(15,:).*y2(3,:)+y2(11,:).*y2(7,:))./r(k,:).^5+4*rho*trapz(x,r(k,:).^2-V)*r(k,:);
    norm_dF(k)=sqrt(trapz(x,dF(k,:).^2));    
    err2=abs((F_prop(k,:)-F_prop(k-1,:))/F_prop(k,:))

    r(k,:)=Proj(r(k-1,:)-sigma*dF(k-1,:)/norm_dF(k-1),a1);
    err1=(max(abs(r(k,:)-r(k-1,:))))
    %err3=sqrt(trapz(x,(r(k,:)-Proj(r(k,:)-dF(k,:),a1,a2)).^2))
end


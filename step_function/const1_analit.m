function [f,g,F_prop,df,dg,F_prop1,F_prop2] = const1_analit(x0,x1,npts,a2)

a=(1/a2)^(1/4)*(cos(pi/8)-1i*sin(pi/8));
I=1i;
% A = [exp(a*x0),exp(-a*x0),exp(a*x0*I),exp(-a*x0*I);...
%     a*exp(a*x0),-a*exp(-a*x0),I*a*exp(a*x0*I),-I*a*exp(-a*x0*I);... % row 2
%     exp(a),exp(-a),-exp(I*a),-exp(-I*a);...
%     exp(a),-exp(-a),-I*exp(I*a),I*exp(-I*a)];

% A = [1,1,1,1;...
%     a,-a,I*a,-I*a;... % row 2
%     exp(a*(1-x0)),exp(-a*(1-x0)),-exp(I*a*(1-x0)),-exp(-I*a*(1-x0));...
%     exp(a*(1-x0)),-exp(-a*(1-x0)),-I*exp(I*a*(1-x0)),I*exp(-I*a*(1-x0))];

A = [1,1,1,1;...
    a,-a,I*a,-I*a;... % row 2
     exp(a),exp(-a),-exp(I*a),-exp(-I*a);...
     exp(a),-exp(-a),-I*exp(I*a),I*exp(-I*a)];

B=[-I*x0/(1-x0);-I;0;0];

C=A\B;



y=C(1)*exp(a*x1)+C(2)*exp(-a*x1)+C(3)*exp(I*a*x1)+C(4)*exp(-I*a*x1);

dy=C(1)*a*exp(a*x1)-C(2)*a*exp(-a*x1)+C(3)*a*I*exp(I*a*x1)-C(4)*a*I*exp(-I*a*x1);

d2y=a^2*(C(1)*exp(a*x1)+C(2)*exp(-a*x1)-C(3)*exp(I*a*x1)-C(4)*exp(-I*a*x1));

d3y=a^3*(C(1)*exp(a*x1)-C(2)*exp(-a*x1)-I*C(3)*exp(I*a*x1)+I*C(4)*exp(-I*a*x1));

f=real(y); g=-imag(y);
df=real(dy); dg=-imag(dy);
F_prop = trapz(x1,df.*g-dg.*f)/2;

F_prop1=a2*(dy(1)*conj(d3y(1))+conj(dy(1))*d3y(1))/4;
F_prop2=-a2*(d2y(1)*conj(d2y(1)))/4;

F_prop = F_prop1+F_prop2;

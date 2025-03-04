function [f,g,F_prop,df,dg,F_prop1,F_prop2] = const_2_analit(x0,npts,a1,a2)


matrices;
x1 = linspace(0,x0,floor(npts*x0));
x2 = linspace(x0,1,npts-floor(npts*x0));
x=[x1,x2];
y1=C(1)*exp(a*x1)+C(2)*exp(-a*x1)+C(3)*exp(I*a*x1)+C(4)*exp(-I*a*x1);
y2=C(5)*exp(b*x2)+C(6)*exp(-b*x2)+C(7)*exp(I*b*x2)+C(8)*exp(-I*b*x2);
y=[y1,y2];
dy1=C(1)*a*exp(a*x1)-C(2)*a*exp(-a*x1)+C(3)*a*I*exp(I*a*x1)-C(4)*a*I*exp(-I*a*x1);
dy2=C(5)*b*exp(b*x2)-C(6)*b*exp(-b*x2)+C(7)*b*I*exp(I*b*x2)-C(8)*b*I*exp(-I*b*x2);

d2y1=a^2*(C(1)*exp(a*x1)+C(2)*exp(-a*x1)-C(3)*exp(I*a*x1)-C(4)*exp(-I*a*x1));
d2y2=b^2*(C(5)*exp(b*x2)+C(6)*exp(-b*x2)-C(7)*exp(I*b*x2)-C(8)*exp(-I*b*x2));

d3y1=a^3*(C(1)*exp(a*x1)-C(2)*exp(-a*x1)-I*C(3)*exp(I*a*x1)+I*C(4)*exp(-I*a*x1));
%d3y2=b^3*(C(5)*exp(b*x2)+C(6)*exp(-b*x2)-C(7)*exp(I*b*x2)-C(8)*exp(-I*b*x2));

dy=[dy1,dy2];
f=real(y); g=-imag(y);
df=real(dy); dg=-imag(dy);
F_prop = trapz(x,df.*g-dg.*f)/2;

F_prop1 = trapz(x1,real(dy1).*(-imag(y1))-(-imag(dy1)).*real(y1))/2;
F_prop2 = trapz(x2,real(dy2).*(-imag(y2))-(-imag(dy2)).*real(y2))/2;

F_prop = (dy1(1)*a2*conj(d3y1(1))+conj(dy1(1))*a2*d3y1(1)-a2*d2y1(1)*conj(d2y1(1))+a2/a1*(a1-a2)*d2y1(end)*conj(d2y1(end)))/4;

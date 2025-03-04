
a=(1/a2)^(1/4)*(cos(pi/8)-1i*sin(pi/8));
b=(1/a1)^(1/4)*(cos(pi/8)-1i*sin(pi/8));
I=1i;
A = [1,1,1,1,0,0,0,0;...
    a,-a,I*a,-I*a,0,0,0,0;... % row 2
    exp(a*x0),exp(-a*x0),exp(a*x0*I),exp(-I*a*x0),...
    -exp(b*x0),-exp(-b*x0),-exp(b*x0*I),-exp(-I*b*x0);... row 3 - cont y
    a*exp(a*x0),-a*exp(-a*x0),a*I*exp(a*x0*I),-a*I*exp(-I*a*x0),...
    -b*exp(b*x0),b*exp(-b*x0),-b*I*exp(b*x0*I),b*I*exp(-I*b*x0);... % row 4 cont y'
    a2*a^2*exp(a*x0),a2*a^2*exp(-a*x0),-a2*a^2*exp(a*x0*I),-a2*a^2*exp(-I*a*x0),...
    -a1*b^2*exp(b*x0),-a1*b^2*exp(-b*x0),a1*b^2*exp(b*x0*I),a1*b^2*exp(-I*b*x0);... % row 5 cont y''
    a2*a^3*exp(a*x0),-a2*a^3*exp(-a*x0),-I*a2*a^3*exp(a*x0*I),I*a2*a^3*exp(-I*a*x0),...
    -a1*b^3*exp(b*x0),a1*b^3*exp(-b*x0),a1*b^3*I*exp(b*x0*I),-a1*b^3*I*exp(-I*b*x0);... % row 6 cont y'''
    0,0,0,0,exp(b),exp(-b),-exp(I*b),-exp(-I*b);...
    0,0,0,0,exp(b),-exp(-b),-I*exp(I*b),I*exp(-I*b)];

B=[0;-I;0;0;0;0;0;0];

C=A\B;


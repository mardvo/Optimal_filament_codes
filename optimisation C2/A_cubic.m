%% function cubic 

function A=A_cubic(a,b,x)

a1=a-b;
a2=-3*a1;
a3=-a2;
a4=b;

A=a1*x.^3+a2*x.^2+a3*x+a4;
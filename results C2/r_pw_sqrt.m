function A=r_pw_sqrt(a1,V,x,x2)

A=zeros([1,length(x)]);

%a=2*(V-a1^2)/x2^2;
b=a1*2/3*x2^(3/2);
a=(-b+sqrt(b^2-x2^2/2*(a1^2-V)))/x2^2*2;
A = a*sqrt(x2-x)+a1;

A(x>=x2)=a1;
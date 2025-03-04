function df = first_der(f,x)
%%%% calculates second order first derivative for uniform arrays
hh=x(2)-x(1); %%% for uniform spacing only
df=zeros(length(x),1); 
df(1)=(-3*f(1)+4*f(2)-f(3))/2/hh;
df(2:length(x)-1)=(f(3:length(x))-f(1:length(x)-2))/2/hh;
df(length(x))=(f(length(x)-2)-4*f(length(x)-1)+3*f(length(x)))/2/hh;

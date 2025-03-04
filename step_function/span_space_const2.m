clear all;
close all;

npts=1000;

ni=100;
a1s=10.^linspace(-3,1,ni);
a2s=10.^linspace(-3,1,ni);
% a1s=0.1;
 a2s=1;
nx0=100;
x0=linspace(0+1/nx0,1-1/nx0,nx0);
F=zeros([ni,ni,nx0]);
for m=1:length(a2s)
    for mm=1:length(a1s)
    a1=a1s(mm);
    a2=a2s(m);
    for i=1:length(x0)
    [~,~,F(mm,m,i),~,~] = const_2_analit(x0(i),npts,a1,a2);
    end
    end
end

save pw_const_a1_1.mat
function dydx = twoode(x,y,A)

f=y(1,:); df=y(2,:); d2f=y(3,:); d3f=y(4,:); g=y(5,:); dg=y(6,:); d2g=y(7,:); d3g=y(8,:);
l=y(9,:); dl=y(10,:); d2l=y(11,:); d3l=y(12,:);
m=y(13,:); dm=y(14,:); d2m=y(15,:); d3m=y(16,:); 

d2f = d2f./A;
d2g = d2g./A;
d4f = -g;
d4g = f;

d2l = d2l./A;
d2m = d2m./A;
d4l = -m-df;
d4m = l+dg;


dydx=[df;d2f;d3f;d4f;dg;d2g;d3g;d4g;dl;d2l;d3l;d4l;dm;d2m;d3m;d4m];




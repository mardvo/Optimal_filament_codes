function res = twobc(ya,yb)
% ya = 0 condizione in x=0 
% yb = 0 condizione in x=1

res = [ya(1);ya(2);yb(3);yb(4);ya(5);ya(6)-1;yb(7);yb(8);ya(9);ya(10);yb(11);yb(12)+yb(1)/2;ya(13);ya(14);yb(15);yb(16)-yb(5)/2];
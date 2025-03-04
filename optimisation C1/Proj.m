function A=Proj(f,a,b)

% projection to a<f<b
for i=1:length(f)
A(i)=max(a,min(f(i),b));
end
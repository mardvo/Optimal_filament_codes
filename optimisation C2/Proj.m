function A=Proj(f,a)

% projection to a<f
for i=1:length(f)
A(i)=max(a,f(i));
end
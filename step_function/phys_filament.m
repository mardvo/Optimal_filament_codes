% figure for realistic forces
close all 
clear all
phys_pars_PPy; % Change here between PDMS and PPy

npts=300;

x = linspace(0,1,npts);

A0=A; 
consts; % solve main problem
const_cong; % solve the conjugate problem
Fin = F0*F_prop;


ni=100; 
a2s = A0;
a1s=a2s*10.^linspace(-1,0,ni); % pick an array 'a' values (or a single value) 

nx0=200;
x0=linspace(5/nx0,1-5/nx0,nx0);
F=zeros([ni,nx0]);
for mm=1:length(a1s)
    a1=a1s(mm);
    a2=1;
    for i=1:length(x0)
    [~,~,F(mm,i),~,~,F1(mm,i),F2(mm,i)] = const_2_analit(x0(i),npts,a1,a2);

    [~,~,F_end(mm,i),~,~] = const1_analit(x0(i),x,npts,a1/(1-x0(i))^4);
    end
end

figure
plot(a1s,max(-real(F*F0)'),'LineWidth',3)
ylabel('F, N')
xlabel('a')
set(gca,'Fontsize',16)

F_opt = max(max(-real(F)'))*F0;

[FF1,ind1] = max(-real(F)');
[FF2,ind2] = max(FF1);

% x0(ind1(ind2))
% Sp = (1-x0(ind1(ind2)))*(a1s(1))^(-1/4)

disp("Force for uniform filament "+num2str(abs(Fin))+" N");
disp("Force for optimal filament "+num2str(abs(F_opt))+" N");
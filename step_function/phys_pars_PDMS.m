% estimate of propulsion for physical swimmer
%% PDMS 
E = 3.86e+6; %Pa

omega = 3.6*2*pi; % Hz

L=1500e-6; % m filament length
r0=5e-6; % m filament thickness

mu = 1.15e-3; % Pa*s

xi_ort = mu*4*pi/(log(L/r0)+1/2);
xi_par = mu*2*pi/(log(L/r0)-1/2);


A0=L^4*xi_ort*omega;

A_dim = E*r0^4*pi;

A = A_dim/A0;

a0 = 0.4; % rad
epsilon = a0; % 

% scale for F
F0 = (xi_ort-xi_par)*epsilon^2*L^2*omega;
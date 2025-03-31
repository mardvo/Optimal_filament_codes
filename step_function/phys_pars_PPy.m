% estimate of propulsion for physical swimmer
%% PPY
E = 100e+6; %Pa

omega = 20*2*pi; % Hz

L=9e-6; % m filament length
r0=100e-9; % m filament thickness

mu = 0.0152; % Pa*s

xi_ort = mu*4*pi/(log(L/r0)+1/2);
xi_par = mu*2*pi/(log(L/r0)-1/2);


A0=L^4*xi_ort*omega;

A_dim = E*r0^4*pi;

A = A_dim/A0;

a0 = 0.4; % rad
epsilon = a0; %

% scale for F
F0 = (xi_ort-xi_par)*epsilon^2*L^2*omega;


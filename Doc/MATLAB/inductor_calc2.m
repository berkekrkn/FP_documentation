clc;
close;
clear;

N = 10;
dc = 2.67e-3;
l = 4.5e-3;
di =  0.4e-3;
delta = 0.01825e-3;

do = di + 2*delta;
eps0 = 8.854e-12;
mu0 = 4*pi*1e-7;
rho_c = 1.72e-8;

p = l/N;
tmp1 = di/p;

kappa = dc/(sqrt(dc^2+l^2));
kappa_prime = l/(sqrt(dc^2+l^2));
[K,E] = ellipke(kappa^2);
K = 4/(3*pi)*1/kappa_prime*((kappa_prime^2)/(kappa^2)*(K-E)+E-kappa);
Ls = mu0*N^2*(pi*(dc/2)^2)/l*K;
ks = 0.4403;
km = 0.2664;
%km = (log(2*pi)-1.5)*(1-0.982889/(N-0.017111))+(-0.16641/N +0.00479/N^2 +0.001772/N^3)*log(N);
L = Ls - mu0*dc/2*N*(ks+km);
fs = (4*rho_c)/(pi*mu0*(0.4e-3)^2)
lw = N*sqrt((pi*dc)^2+(l/N)^2) + 2*3.39e-3;
Rdc = rho_c*lw/(pi*(di/2)^2);

Cs = 4*eps0*l/pi*(0.71*dc/l+1+2.4*(dc/l)^1.5);
fprintf("L: %.3f nH\nCs: %.3f pF\nRdc: %.3f mOhm\nfSR: %.3f GHz\n", ...
            L*1e9, Cs*1e12, Rdc*1e3, 1/(2*pi*sqrt(L*Cs))*1e-9);
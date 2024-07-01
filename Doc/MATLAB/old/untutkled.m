a = 0.161;
d = 0.021;

MM = (8*(a*log((a+sqrt(a^2+d^2))/(a+sqrt(2*a^2+d^2))*sqrt(a^2+d^2)/d))+ ...
8*(sqrt(2*a^2+d^2)-2*sqrt(a^2+d^2)+d));

close all
clear
clc

%% Inductance

D = 1.59e-3;
d_i = 0.16e-3;
delta = 0.01e-3;
xi_r = -0.005e-3;
xi_z = 0.004e-3;
N = 12;
d_c = 5e-3;

R_i = D/2;
d_o = d_i + 2*delta;
N_V = 6;
N_L = 2;
a = d_o/2;
i = 1:N_L;
j = 1:N_V;
r = R_i+(i-1)*(2*a-xi_r)+a-xi_r;
z = (j-1)*(2*a+xi_z)+a;


% function L = self_ind(b, a)
%     rho = a*1e2;
%     a = b*1e2;
%     L = 8*a*(log(a/rho)+rho/a-0.524)*1e-9;
% end
% function M = mutual_ind(b,c,z)
% mu_0 = 4*pi*1e-7;
%     a = b/2;
%     c = c/2;
%     M = 2*mu_0/pi*(sqrt(2*(a+c)^2+z^2)+sqrt(2*(a-c)^2+z^2) ...
%     - 2*sqrt(2*a^2+2*c^2+z^2)-(a+c)*atanh((a+c)/sqrt(2*(a+c)^2+z^2))...
%     -(a-c)*atanh((a-c)/(sqrt(2*(a-c)^2+z^2))) ...
%     +(a+c)*atanh((a+c)/sqrt(2*a^2+2*c^2+z^2)) ...
%     +(a-c)*atanh((a-c)/sqrt(2*a^2+2*c^2+z^2)));
% end

function M = mutual_ind(r1,r2,x)
    mu_0 = 4*pi*1e-7;
    k = (2*sqrt(r1*r2)/(sqrt((r1+r2)^2+x^2)));
    %x1 = sqrt((r1+r2)^2+x^2);
    %x2 = sqrt((r1-r2)^2+x^2);
    %k = (x1-x2)/(x1+x2);
    [K,E] = ellipke(k^2);
    M = 2*mu_0*sqrt(r1*r2)/k*((1-1/2*k^2)*K-E);
    %M = 2*mu_0*sqrt(r1*r2)/sqrt(k)*(K-E);
end

function L = self_ind(r, rw)
    r1 = r; r2 = r;
    x = rw*exp(-1/4);
    mu_0 = 4*pi*1e-7;
    k = (2*sqrt(r1*r2)/(sqrt((r1+r2)^2+x^2)));
    [K,E] = ellipke(k^2);
    L = -mu_0*sqrt(r1*r2)*((k-2/k)*K+2/k*E);
end


l0 = 0;
cntr1 = 0;
lval = [];
for i = 1:N_L
    for j = 1:N_V
        li = self_ind(r(i), a);
        l0 = l0 + li;
        lval = [lval, li];
        cntr1 = cntr1+1;
    end
end

m = 0;
cntr1 = 0;
cntr2 = 0;
M_vals = [];
for i = 1:N_L 
    for j = 1:N_V % (i,j) of chosen winding
        for k = 1:N_L
            for l = 1:N_V % (k,l) of compared winding
                if i ~= k || j ~= l
                    m_ijkl = mutual_ind(r(i), r(k), abs(z(j)-z(l)));
                    M_vals = [M_vals, m_ijkl];
                    m = m + m_ijkl;
                    cntr2 = cntr2+1;
                else 
                    M_vals = [M_vals, 0];
                end
                cntr1 = cntr1+1;
            end
        end
    end
end

M_vals = M_vals*1e9;
L = l0 + m;
disp(m)

M_real = [1.79752, 1.0275, 0.6745, 0.4666, ...
        2.1197, 1.656, 1.145, 0.8003, 0.4339, ...
        1.4583, 1.33132, 1.0768, 0.644089, 0.4989, ...
        1.58268, 1.0978, 0.79834, 0.59805];
M = reshape(M_vals, 15, 15);
err = abs(M_real - M);
disp(err);
disp(sum(err));
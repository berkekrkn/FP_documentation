close all
clear
clc

%% Inductance

D = 1.4e-3;
d_i = 0.18e-3;
delta = 0.01e-3;
xi_r = -0.005e-3;
xi_z = 0.01e-3;
N = 15;
d_c = 5e-3;

R_i = D/2;
d_o = d_i + 2*delta;
N_V = 5;
N_L = 3;
a = d_o/2;
i = 1:N_L;
j = 1:N_V;

r = R_i+(i-1)*(2*a-xi_r)+a-xi_r;
z = (j-1)*(2*a+xi_z)+a;


function L = self_ind(b, a)
    rho = a*1e2;
    a = b*1e2;
    L = 8*a*(log(a/rho)+rho/a-0.524)*1e-9;
end

function M = mutual_ind(a,c,z1, z2)
    z = abs(z1-z2);
    a = a/2;
    c = c/2;
    mu0 = 4*pi*1e-7;
    alpha = 2*sqrt(a*c/((a+c)^2+z^2));
    [K,E] = ellipke(alpha^2);
    M = 1.1*(2*mu0*sqrt(a*c)/alpha*((1-alpha^2/2)*K-E));
end

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

% function M = mutual_ind_corrected(r1,r2,z1,z2,di)
%     R1 = r1;
%     R2 = r2;
%     Z1 = z1;
%     Z2 = z2;
%     Square_side = di/sqrt(2);
%     N = 15;
%     qt = linspace(-Square_side/2, Square_side/2, N);
%     pt = linspace(-Square_side/2, Square_side/2, N);
%     qr = linspace(-Square_side/2, Square_side/2, N);
%     pr = linspace(-Square_side/2, Square_side/2, N);
%     Qt = length(qt);
%     Pt = length(pt);
%     Qr = length(qr);
%     Pr = length(pr);
%     Ms = 0;
%     for i = 1:Qt
%         for j = 1:Pt
%             for k = 1:Qr
%                 for l = 1:Pr
%                     r1 = R1+pr(l);
%                     z1 = Z1+qr(k);
%                     r2 = R2+pt(j);
%                     z2 = Z2+qt(i);
%                     Mij = mutual_ind(2*r1,2*r2, abs(z1-z2));
%                     Ms = Ms + Mij;
%                 end
%             end
%         end
%     end
%     M = Ms/(Qt*Pt*Qr*Pr);
% end

% function M = mutual_ind(r1,r2,x)
%     mu_0 = 4*pi*1e-7;
%     k = (2*sqrt(r1*r2)/(sqrt((r1+r2)^2+x^2)));
%     %x1 = sqrt((r1+r2)^2+x^2);
%     %x2 = sqrt((r1-r2)^2+x^2);
%     %k = (x1-x2)/(x1+x2);
%     [K,E] = ellipke(k^2);
%     M = 2*mu_0*sqrt(r1*r2)/k*((1-1/2*k^2)*K-E);
%     %M = 2*mu_0*sqrt(r1*r2)/sqrt(k)*(K-E);
% end
% 
% function L = self_ind(r, rw)
%     r1 = r; r2 = r;
%     x = rw*exp(-1/4);
%     mu_0 = 4*pi*1e-7;
%     k = (2*sqrt(r1*r2)/(sqrt((r1+r2)^2+x^2)));
%     [K,E] = ellipke(k^2);
%     L = -mu_0*sqrt(r1*r2)*((k-2/k)*K+2/k*E);
% end


l0 = 0;
cntr1 = 0;
lval = [];
for i = 1:N_L
    for j = 1:N_V
        li = self_ind(2*r(i), a);
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
                    m_ijkl = mutual_ind(2*r(i), 2*r(k), z(j), z(l));
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

M_real =   [2.7130, 1.7975, 1.0275, 0.6745, 0.4666, ...
            2.1197, 1.6560, 1.1450, 0.8083, 0.5858, ...
            1.4582, 1.3313, 1.0769, 0.8357, 0.6441; ...
            1.7975, 2.7130, 1.7975, 1.0275, 0.6745, ...
            1.6560, 2.1197, 1.6560, 1.1450, 0.8083, ...
            1.3313, 1.4582, 1.3313, 1.0769, 0.8357; ...
            1.0275, 1.7975, 2.7130, 1.7975, 1.0275, ...
            1.1450, 1.6560, 2.1197, 1.6560, 1.1450, ...
            1.0769, 1.3313, 1.4582, 1.3313, 1.0769; ...   
            0.6745, 1.0275, 1.7975, 2.7130, 1.7975, ...
            0.8083, 1.1450, 1.6560, 2.1197, 1.6560, ...
            0.8357, 1.0769, 1.3313, 1.4582, 1.3313; ...
            0.4666, 0.6745, 1.0275, 1.7975, 2.7130, ...
            0.5858, 0.8083, 1.1450, 1.6560, 2.1197, ...
            0.6441, 0.8357, 1.0769, 1.3313, 1.4582; ...
            2.1197, 1.6560, 1.1450, 0.8083, 0.5858, ...
            3.7640, 2.5956, 1.5827, 1.0978, 0.7983, ...
            2.9580, 2.3644, 1.6987, 1.2483, 0.9412; ...
            1.6560, 2.1197, 1.6560, 1.1450, 0.8083, ...
            2.5956, 3.7640, 2.5956, 1.5827, 1.0978, ...
            2.3644, 2.9580, 2.3644, 1.6987, 1.2483; ...
            1.1450, 1.6560, 2.1197, 1.6560, 1.1450, ...
            1.5827, 2.5956, 3.7640, 2.5956, 1.5827, ...
            1.6987, 2.3644, 2.9580, 2.3644, 1.6987; ...
            0.8083, 1.1450, 1.6560, 2.1197, 1.6560, ...
            1.0978, 1.5827, 2.5956, 3.7640, 2.5956, ...
            1.2483, 1.6987, 2.3644, 2.9580, 2.3644; ...
            0.5858, 0.8083, 1.1450, 1.6560, 2.1197, ...
            0.7983, 1.0978, 1.5827, 2.5956, 3.7640, ...
            0.9412, 1.2483, 1.6987, 2.3644, 2.9580; ...
            1.4582, 1.3313, 1.0769, 0.8357, 0.6441, ...
            2.9580, 2.3644, 1.6987, 1.2483, 0.9412, ...
            4.8580, 3.4559, 2.2039, 1.5866, 1.1953; ...
            1.3313, 1.4582, 1.3313, 1.0769, 0.8357, ...
            2.3644, 2.9580, 2.3644, 1.6987, 1.2483, ...
            3.4559, 4.8580, 3.4559, 2.2039, 1.5866; ...
            1.0769, 1.3313, 1.4582, 1.3313, 1.0769, ...
            1.6987, 2.3644, 2.9580, 2.3644, 1.6987, ...
            2.2039, 3.4559, 4.8580, 3.4559, 2.2039; ...
            0.8357, 1.0769, 1.3313, 1.4582, 1.3313, ...
            1.2483, 1.6987, 2.3644, 2.9580, 2.3644, ...
            1.5866, 2.2039, 3.4559, 4.8580, 3.4559; ...
            0.6441, 0.8357, 1.0769, 1.3313, 1.4582, ...
            0.9412, 1.2483, 1.6987, 2.3644, 2.9580, ...
            1.1953, 1.5866, 2.2039, 3.4559, 4.8580];
M = reshape(M_vals, 15, 15);
M = M + diag(lval);
err = abs(M_real - M);
disp(err);
disp(sum(err));
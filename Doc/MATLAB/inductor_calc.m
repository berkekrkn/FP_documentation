close all
clear
clc

%Input parameters
D = 1.4e-3; % Rod side length (m)
d_i = 0.18e-3; % Wire inner radius (m)
delta = 0.01e-3; % Insulation thickness (m)
xi_r = -0.005e-3; % Negative spacing between adjacent layers (m)
xi_z = 0.009e-3; % Spacing between adjacent windings in the same layer (m)
N_V = 5; % Number of windings per layer
N_L = 3; % Number of layers

% D = 1.4e-3; % Rod side length (m)
% d_i = 0.16e-3; % Wire inner radius (m)
% delta = 0.01e-3; % Insulation thickness (m)
% xi_r = -0.005e-3; % Negative spacing between adjacent layers (m)
% xi_z = 0.005e-3; % Spacing between adjacent windings in the same layer (m)
% N_V = 6; % Number of windings per layer
% N_L = 2; % Number of layers

% Intermediate parameters
N = N_V*N_L; % Total number of windings
R_i = D/2; % Self explanatory
d_o = d_i + 2*delta; % Wire outer diameter (m)
d_c = (N_V-1)*xi_z + N_V*d_o; % Coil length (m)
a = d_o/2; % Wire outer radius (m)

% Coil model
i = 1:N_L;
j = 1:N_V;

% Radius of the winding center points from the cylindrical axis
r = R_i+(i-1)*(2*a-xi_r)+a-xi_r; 
% Distance of middle points of windings from the zero point
z = (j-1)*(2*a+xi_z)+a;

function L = self_ind(a, rho)
    % Calculates self inductance of a square loop
    % a: Loop side length
    % rho: Loop wire radius
    mu_0 = 4*pi*1e-7;
    L = mu_0*(2*a/pi)*(log(a/rho) - 0.77401);
end

function M = mutual_ind(b,c,z)
    % Calculates mutual inductance of two square loops with side length
    % b and c, which are z apart
    % b: Side length of loop #1 (m)
    % c: Side length of loop #2 (m)
    % z: Distance between two loops (m)
    mu_0 = 4*pi*1e-7;
    a = b/2;
    c = c/2;
    M = 2*mu_0/pi*(sqrt(2*(a+c)^2+z^2)+sqrt(2*(a-c)^2+z^2) ...
    - 2*sqrt(2*a^2+2*c^2+z^2)-(a+c)*atanh((a+c)/sqrt(2*(a+c)^2+z^2))...
    -(a-c)*atanh((a-c)/(sqrt(2*(a-c)^2+z^2))) ...
    +(a+c)*atanh((a+c)/sqrt(2*a^2+2*c^2+z^2)) ...
    +(a-c)*atanh((a-c)/sqrt(2*a^2+2*c^2+z^2)));
end

% Sum the self-inductance of each windings
l0 = 0;
cntr1 = 0;
lval = [];
for i = 1:N_L
    for j = 1:N_V
        li = self_ind(2*r(i), d_i/2);
        l0 = l0 + li;
        lval = [lval, li];
        cntr1 = cntr1+1;
    end
end

% Sum the muutual inductance of each winding pairs
m = 0;
cntr1 = 0;
cntr2 = 0;
M_vals = [];
for i = 1:N_L 
    for j = 1:N_V % (i,j) is the chosen winding
        for k = 1:N_L
            for l = 1:N_V % (k,l) is the winding that is being compared to (i,j)
                if i ~= k || j ~= l
                    m_ijkl = mutual_ind(2*r(i), 2*r(k), abs(z(j)-z(l)));
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

%% Resistance calculation
rho_c = 1.72e-8; % Copper resistance
p = d_o+xi_z; % Coil pitch
lw = sum(4*N_V*sqrt((2*r).^2+(p/4)^2));
Rdc = rho_c/(pi*(d_i/2)^2)*lw;

% Capacitance of multilayer coil calculation: DEPRECATED
% D = 1.98e-3;
% R_i = D/2;
% eps_0 = 8.854e-12;
% eps_r = 1;
% d_tilde = 1.26*d_o-1.15*d_i+xi_z;
% Cs = 8*pi*eps_0*eps_r*d_c*(N_L-1)/(6*d_tilde*N_L^2)*(2*R_i+d_o-N_L*xi_r);
% fSR = 1/(2*pi*sqrt(1.5e-6*Cs));

%% Error Analysis for 744045002 (comment out if other coil dimenasions used)
% Real values from CST simulations (For N_V = 5, N_L = 3 coil
% M_real =   [2.7130, 1.7975, 1.0275, 0.6745, 0.4666, ...
%             2.1197, 1.6560, 1.1450, 0.8083, 0.5858, ...
%             1.4582, 1.3313, 1.0769, 0.8357, 0.6441; ...
%             1.7975, 2.7130, 1.7975, 1.0275, 0.6745, ...
%             1.6560, 2.1197, 1.6560, 1.1450, 0.8083, ...
%             1.3313, 1.4582, 1.3313, 1.0769, 0.8357; ...
%             1.0275, 1.7975, 2.7130, 1.7975, 1.0275, ...
%             1.1450, 1.6560, 2.1197, 1.6560, 1.1450, ...
%             1.0769, 1.3313, 1.4582, 1.3313, 1.0769; ...   
%             0.6745, 1.0275, 1.7975, 2.7130, 1.7975, ...
%             0.8083, 1.1450, 1.6560, 2.1197, 1.6560, ...
%             0.8357, 1.0769, 1.3313, 1.4582, 1.3313; ...
%             0.4666, 0.6745, 1.0275, 1.7975, 2.7130, ...
%             0.5858, 0.8083, 1.1450, 1.6560, 2.1197, ...
%             0.6441, 0.8357, 1.0769, 1.3313, 1.4582; ...
%             2.1197, 1.6560, 1.1450, 0.8083, 0.5858, ...
%             3.7640, 2.5956, 1.5827, 1.0978, 0.7983, ...
%             2.9580, 2.3644, 1.6987, 1.2483, 0.9412; ...
%             1.6560, 2.1197, 1.6560, 1.1450, 0.8083, ...
%             2.5956, 3.7640, 2.5956, 1.5827, 1.0978, ...
%             2.3644, 2.9580, 2.3644, 1.6987, 1.2483; ...
%             1.1450, 1.6560, 2.1197, 1.6560, 1.1450, ...
%             1.5827, 2.5956, 3.7640, 2.5956, 1.5827, ...
%             1.6987, 2.3644, 2.9580, 2.3644, 1.6987; ...
%             0.8083, 1.1450, 1.6560, 2.1197, 1.6560, ...
%             1.0978, 1.5827, 2.5956, 3.7640, 2.5956, ...
%             1.2483, 1.6987, 2.3644, 2.9580, 2.3644; ...
%             0.5858, 0.8083, 1.1450, 1.6560, 2.1197, ...
%             0.7983, 1.0978, 1.5827, 2.5956, 3.7640, ...
%             0.9412, 1.2483, 1.6987, 2.3644, 2.9580; ...
%             1.4582, 1.3313, 1.0769, 0.8357, 0.6441, ...
%             2.9580, 2.3644, 1.6987, 1.2483, 0.9412, ...
%             4.8580, 3.4559, 2.2039, 1.5866, 1.1953; ...
%             1.3313, 1.4582, 1.3313, 1.0769, 0.8357, ...
%             2.3644, 2.9580, 2.3644, 1.6987, 1.2483, ...
%             3.4559, 4.8580, 3.4559, 2.2039, 1.5866; ...
%             1.0769, 1.3313, 1.4582, 1.3313, 1.0769, ...
%             1.6987, 2.3644, 2.9580, 2.3644, 1.6987, ...
%             2.2039, 3.4559, 4.8580, 3.4559, 2.2039; ...
%             0.8357, 1.0769, 1.3313, 1.4582, 1.3313, ...
%             1.2483, 1.6987, 2.3644, 2.9580, 2.3644, ...
%             1.5866, 2.2039, 3.4559, 4.8580, 3.4559; ...
%             0.6441, 0.8357, 1.0769, 1.3313, 1.4582, ...
%             0.9412, 1.2483, 1.6987, 2.3644, 2.9580, ...
%             1.1953, 1.5866, 2.2039, 3.4559, 4.8580];
% %Error analysis between CST and analytical equations
% M = reshape(M_vals, N, N);
% M = M + diag(lval)*1e9;
% errAbs = abs(M_real - M);
% errRel = abs(M_real - M)./M_real*100;
% sum(M_real, "all")
% sum(M, "all")
% sum(errAbs, "all")
% mean(errRel, "all")
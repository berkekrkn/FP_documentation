close all;
clc;
clear;

data = readmatrix("data/15.csv");
init_f = linspace(0, 1e3, 1001)';
omega = 2*pi*[init_f; data(:,1)];
init_mureal = 1582*ones(1,length(init_f))';
init_muim = 1i*zeros(1,length(init_f))';
mu_complex = complex([init_mureal; data(:,2)], ...
                    [init_muim; -data(:,3)]);
mu_complex2 = complex(data(:,2), -data(:,3));
f2 =  data(:,1);
f = logspace(0,10,10001);
funcdw = @(p,omega) (p(1)*p(2)^2)./(p(2)^2-omega.^2+1i*omega*p(3));
funcsp = @(p,omega) (p(1)*p(2)*(p(2)+1i*omega*p(3)))./((p(2)+1i*omega*p(3)).^2-omega.^2);
func = @(p,omega) 1 + funcdw([p(1), p(2), p(3)], omega) + funcsp([p(4),p(5),p(6)], omega);
objFunc = @(p) sum(abs([real(func(p, omega)) - real(mu_complex); ...
                            imag(func(p, omega)) - imag(mu_complex)]));
%dw: p(1): X_dw, p(2): omega_dw, p(3): beta
lbdw = [0, 1e8, 1e9];
ubdw = [1500, 1e10, 1e12];
%sp: p(4): X_sp, p(5): omega_sp, p(6): alpha
lbsp = [0, 1e8, 0];
ubsp = [1500, 1e11, 500];
lb = [lbdw, lbsp];
ub = [ubdw, ubsp];

p0 = [117, 1e9, 6.5e9, 1450, 7e9, 330];
%options = optimoptions('simulannealbnd', 'Display', 'iter');
options = optimoptions('particleswarm', 'Display','iter', ...
    'MaxStallIterations',50, 'FunctionTolerance', 1e-6);
%options = optimoptions( 'ga', 'Display', 'iter',...
             %'MaxStallGenerations', 100);
bestval = inf;
best = zeros(1,6);
for i = 1:1
    %[bestParams, fval, exitflag, output] = simulannealbnd(objFunc, p0, lb, ub, options);
    [bestParams, fval, exitflag, output] = particleswarm(objFunc, 6, lb, ub, options);
    %[bestParams, fval, exitflag, output] = ga(objFunc, 6, [], [], [],[], lb, ub, [], options);
    if fval < bestval
        best = bestParams;
        bestval = fval;
    end
end
p02 = best;
options = optimoptions('fmincon', 'Display', 'iter', ...
    'Algorithm', 'interior-point');
best = fmincon(objFunc,best,[],[],[],[],lb,ub,[],options);
%%
figure;
plot(f2/1e6, real(mu_complex2), '.', 'LineWidth', 1.5);
set(gca, 'XScale', 'Log');
hold on; grid on;
plot(f/1e6, real(func(best, 2*pi*f)), '-', 'LineWidth', 1.5);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'Log');
xlabel('Frequnecy, f in MHz');
ylabel('Real part of $\underline{\mu_r}$', 'Interpreter', 'latex');
legend('Measurement', 'Curve fit', 'FontSize', 14, 'location', 'best');
xticks([1e-6; 1e-4; 1e-2; 1e0; 1e2; 1e4]);

figure;
plot(f2/1e6, -imag(mu_complex2), '.', 'LineWidth', 1.5, 'Color', "#7E2F8E");
set(gca, 'XScale', 'Log');
hold on; grid on;
plot(f/1e6, -imag(func(best, 2*pi*f)), '-', 'LineWidth', 1.5, 'Color', "#77AC30");
legend('Measurement', 'Curve fit', 'FontSize', 14, 'location', 'best');
xticks([1e-6; 1e-4; 1e-2; 1e0; 1e2; 1e4]);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'Log');
xlabel('Frequnecy, f in MHz');
ylabel('Imag. part of $\underline{\mu_r}$', 'Interpreter', 'latex');

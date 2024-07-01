close all
clear
clc
%% Demagnetization factor
m = 2:0.01:40;
demag = 0.37*m.^(-1.44)*10^(5e-2);
figure;
grid on;
hold on;
semilogx(m, demag, 'LineWidth', 1.5);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'log', 'YScale', 'log');
tickLocx = [2, 3, 4, 5, 10, 20, 30, 40];
set(gca, 'XLim',[2 40], 'XTick', tickLocx, 'XTickLabel',tickLocx);
xlabel('Length-to-diameter ratio');
ylabel('Demagnetization factor, D');
ylim([1e-4, 1]);
%% Rod permeability vs. length-to-diameter ratio and physical permeability
figure;
grid on;
hold on;
mur = [10, 25, 50, 100, 250, 500, 1000];
murod = zeros(length(m),length(mur));
for i = 1:length(mur)
    murod(:,i) = mur(i)./(1+demag*(mur(i)-1));
end
for i = 1:length(mur)
    plot(m, murod(:,i), 'LineWidth', 1.5, 'Color', "#0072BD");
    if i == 1
        text(m(end), murod(end,i), "\mu_r = " + num2str(mur(i)), ...
        'FontSize',16, 'FontName', 'Times New Roman');
    else
        text(m(end), murod(end,i), num2str(mur(i)), ...
        'FontSize', 16, 'FontName', 'Times New Roman');
    end
end
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
xlabel('Length-to-diameter ratio');
ylabel('Eff. rel. permeability, \mu_{r,eff}');
tickLocx = [2, 3, 4, 5, 10, 20, 30, 40];
tickLocy = [10, 25, 50, 100, 250, 500];
set(gca, 'XLim',[2, 40], 'XTick', tickLocx, 'XTickLabel',tickLocx);
set(gca, 'YLim',[1, 500], 'YTick', tickLocy, 'YTickLabel',tickLocy);
%% Fair-Rite Temp vs mu_r
tmp = readmatrix("data/temp.csv");
mur = tmp(:,2);
temp = tmp(:,1);
smoothdata(temp, "gaussian", 100);
tempi = temp(1):0.01:temp(end);
muri = interp1(temp, mur, tempi,"spline");
figure;
plot(temp, mur, 'LineWidth', 1.5);
ylim([0, 2500]);
xlabel('Temperature in $^\circ \mathrm{C}$', 'Interpreter','latex');
ylabel('Relative permeability, \mu_r');
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
%% 7440450015 inductance vs. current
L015 = readmatrix("data/7440450015.csv");
I = L015(1,1):0.01:L015(end,1);
L015_interp = interp1(L015(:,1),L015(:,2),I,"spline"); 
figure
grid on;
hold on;
plot(I, L015_interp,'LineWidth', 1.5);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
xlabel('Current, I in A');
ylabel('Inductance, L in µH');
ylim([0, 1.6]);
xlim([0, 10]);
%% Magnetization curves
par = [1.6, 0.7, 5, 0, -1, 0.01];
syms sig(x) dsig(x) sig2(x)
sig(x) = par(1)/(par(2)+exp(-par(3)*x+par(4)))+par(5)+par(6)*x;
sig2(x) = 1/(1+exp(-x));
dsig(x) = diff(sig(x), x);
xi = logspace(-5, 2, 1001);
xi2 = -7.5:0.01:7.5;
xi3 = 0:0.01:7.5;
y = double(sig(xi))./xi;
y(1,1:551) = 2;
figure;
hold on;
p = plot(xi2, double(sig2(xi2-2)), 'LineWidth', 1.5, 'Color', "#0072BD");
plot(xi2, double(sig2(xi2+2)), 'LineWidth', 1.5, 'Color', "#0072BD");
plot(xi3, double(sig2(xi3)), '--', 'LineWidth', 1.5, 'Color', "#0072BD");
p.Parent.XAxisLocation = "origin";
p.Parent.YAxisLocation ="origin";
xticks(0);
yticks(0);
set(gca,'XColor', 'none','YColor','none')
figure;
hold on;
%plot(xi, sig(xi), 'Linewidth', 1.5);
plot(xi, y, 'Linewidth', 1.5);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'log');
%xlim([1e-5, 1e2]);
xlabel("Magnetic field strength, H in A/m");
ylabel("Rel. permeability, \mu_{r}");
yticks([0, 2.51]);
yticklabels({'1', '\mu_{r,i}'});
xticks([1e-5, 1]);
xticklabels({'10^{-5}H_{sat}', '', ...
           '', 'H_{sat}'});
%% Fair-rite permeability vs. frequency 
mu = readmatrix("data/15.csv");
freq = mu(:,1);
mu_real = mu(:,2);
mu_im = mu(:,3);
figure;
grid on;
hold on;
plot(freq/1e6, mu_real, 'LineWidth', 2);
plot(freq/1e6, mu_im, 'LineWidth', 2);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'log');
xticks([1e-2, 1e-1, 1e0, 1e1, 1e2]);
%ylim([-10 200]);
%xlim([1 1000]);
%curtick = get(gca, 'XTick');
%set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
legend('Real part: \mu_r''', 'Imag. part: \mu_r''''');
xlabel('Frequency, f in MHz');
ylabel('Complex rel. permeability, $\underline{\mu_r}$', 'Interpreter', 'latex');
%% Debye curves
mu_i = 125;
mu_inf = -5;
tau = 4e-9;
freq = logspace(0, 9, 10000);
mu_deb = (mu_i-mu_inf)./(1+1i*(2*pi*tau*freq).^2)+mu_inf;
figure;
hold on;
grid on;
plot(freq/1e6, real(mu_deb), 'LineWidth', 1.5);
plot(freq/1e6, -imag(mu_deb), 'LineWidth', 1.5);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'log');
ylim([-10 200]);
xlim([1 1000]);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
legend('Real part: \mu_r''', 'Imag. part: \mu_r''''');
xlabel('Frequency, f in MHz');
ylabel('Complex rel. permeability, $\underline{\mu_r}$', 'Interpreter', 'latex');
%% 7440450015 impedance vs. freq comparison 
magZ_0015_cst = readmatrix("data/magZ_0015_cst.txt");
magZ_0015_lt = readmatrix("data/magZ_0015_lt.txt");
freq_cst = magZ_0015_cst(:,1);
magZ_0015_imp_cst = magZ_0015_cst(:,2);
freq_lt = magZ_0015_lt(:,1)*1e-6;
magZ_0015_imp_lt = magZ_0015_lt(:,2);
magZ_0015_imp_lt = 10.^(magZ_0015_imp_lt./20);
figure;
grid on;
hold on;
plot(freq_cst, magZ_0015_imp_cst, 'LineWidth', 1.5);
plot(freq_lt, magZ_0015_imp_lt, 'LineWidth', 1.5, 'LineStyle','--');
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'log');
ylim([0 5000]);
xlim([1e-3 1e3]);
xticks([1e-3, 1e-2, 1e-1, 1, 10, 100, 1000])
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
legend('CST Simulation', 'Manufacturer Data', 'location', 'best');
xlabel('Frequency, f in MHz');
ylabel('Magnitude of impedance in $\Omega$', 'Interpreter', 'latex');
%% 744045002 impedance vs. freq comparison  
magZ_002_cst = readmatrix("data/magZ_002_cst.txt");
magZ_002_lt = readmatrix("data/magZ_002_lt.txt");
freq_cst = magZ_002_cst(:,1);
magZ_002_imp_cst = magZ_002_cst(:,2);
freq_lt = magZ_002_lt(:,1)*1e-6;
magZ_002_imp_lt = magZ_002_lt(:,2);
magZ_002_imp_lt = 10.^(magZ_002_imp_lt./20);
figure;
grid on;
hold on;
plot(freq_cst, magZ_002_imp_cst, 'LineWidth', 1.5);
plot(freq_lt, magZ_002_imp_lt, 'LineWidth', 1.5, 'LineStyle','--');
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'log');
ylim([0 10000]);
xlim([1e-3 1e3]);
xticks([1e-3, 1e-2, 1e-1, 1, 10, 100, 1000])
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
legend('CST Simulation', 'Manufacturer Data', 'location', 'best');
xlabel('Frequency, f in MHz');
ylabel('Magnitude of impedance in $\Omega$', 'Interpreter', 'latex');
%% SMD Resistor impedance vs. freq
magZ_res = readmatrix("data/res.txt");
freq_res = magZ_res(:,1);
magZ_res_imp = magZ_res(:,2);
figure;
grid on;
hold on;
plot(freq_res, magZ_res_imp, 'LineWidth', 1.5);
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'log');
ylim([0 100]);
xlim([1e-2 5e3]);
xticks([1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 5000])
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));
xlabel('Frequency, f in MHz');
ylabel('Magnitude of impedance in $\Omega$', 'Interpreter', 'latex');











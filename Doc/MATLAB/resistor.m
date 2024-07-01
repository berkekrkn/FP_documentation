close all;
clear;
clc;


%% Functions
function [m, mdb] = read_mag(matrix)
    x = -200:1:200;
    y = -200:1:200;
    [X,~] = meshgrid(x,y);
    Mx = complex(matrix(:,4), matrix(:,5));
    My = complex(matrix(:,6), matrix(:,7));
    Mz = complex(matrix(:,8), matrix(:,9));
    M = [Mx, My, Mz];
    Mmag = reshape(vecnorm(M,2,2),size(X)).';
    m = Mmag; 
    mdb = 20*log10(Mmag);
end

function plot_field(field, type, freq)
    x = -200:1:200;
    y = -200:1:200;
    [X,Y] = meshgrid(x,y);
    ax = gca;
    h = surf(X,Y,field,'Parent',ax);
    set(gca,'FontSize', 12);
    set(gca,'FontName', "Times New Roman");
    set(h, 'edgecolor','none');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    view(ax,[0,90]);
    colormap(hot);
    c = colorbar('orientation','horizontal','Location','southoutside');
    c.FontSize = 10;
    if type == 'e'
        clim([-20,20]);
        c.Ticks = [-20,20];
    elseif type == 'h'
        clim([-60,0]);
        c.Ticks = [-60,0];
    end
    if 0 < freq && freq < 1000
        title(num2str(freq)+" Hz", "FontName", "Times New Roman", ...
            "FontSize", 12);
    elseif 1000 <= freq && freq < 1e6
        title(num2str(freq/1e3)+" kHz", "FontName", "Times New Roman", ...
            "FontSize", 12);
    elseif 1e6 <= freq && freq < 1e9
        title(num2str(freq/1e6)+" MHz", "FontName", "Times New Roman", ...
            "FontSize", 12);
    elseif 1e9 <= freq
        title(num2str(freq/1e9)+" GHz", "FontName", "Times New Roman", ...
            "FontSize", 12);
    end
end


freqs = ["1e-6", "0.0001", "0.01", "1", "10", "100", "500", "1000"];
allfreqs = ["1e-6", "1e-5", "0.0001", "0.001", "0.01", "0.1", "1", ...
            "10", "30", "50", "100", "250", "500", "750", "1000", ...
            "1500", "2000", "2500", "3000"];
N = length(freqs);
freqstr = "(f=" + freqs + ")";
%% Resistor
 pe = "data/fields/resistor/e/";
 e = dir(fullfile(pe,'*.txt'));  
 fig = figure("Position", [575, 400, 800, 200]);
for i = 1:numel(e)
    f = fullfile(pe,e(i).name);
    if contains(f, freqstr)
        data = readmatrix(f);
        s1 = extractAfter(f,'=');
        s2 = extractBefore(s1,')');
        freq = str2double(s2);
        s = num2str(freq);
        if s(s=='e') 
            s(s=='0') = []; 
        end
        idx = find(freqs == s);
        subplot(1,N,idx);
        [~, mdb] = read_mag(data);
        plot_field(mdb, 'e', freq*1e6);
    end
end
sgtitle('Electric Field Magnitude in dBV/m','FontName','Times New Roman', ...
        'FontSize', 14);
figure("Position", [575, 400, 800, 200]);
ph = "data/fields/resistor/h/";
h = dir(fullfile(ph,'*.txt'));
for i = 1:numel(h)
    f = fullfile(ph,h(i).name);
    if contains(f, freqstr)
        data = readmatrix(f);
        s1 = extractAfter(f,'=');
        s2 = extractBefore(s1,')');
        freq = str2double(s2);
        s = num2str(freq);
        if s(s=='e') 
            s(s=='0') = []; 
        end
        idx = find(freqs == s);
        subplot(1,8,idx);
        [~, mdb] = read_mag(data);
        plot_field(mdb, 'h', freq*1e6);
    end
end
sgtitle('Magnetic Field Magnitude in dBA/m','FontName','Times New Roman', ...
        'FontSize', 14);

avgE = zeros(1,numel(e));
avgH = zeros(1,numel(h));
avgZw = zeros(1,numel(h));
freq = str2double(allfreqs);
for i = 1:numel(e)
    f = fullfile(pe,e(i).name);
    data = readmatrix(f);
    [mage, ~] = read_mag(data);
    s1 = extractAfter(f,'=');
    s2 = extractBefore(s1,')');
    avgE(allfreqs == s2) = mean(mage, "all");
    f = fullfile(ph,h(i).name);
    data = readmatrix(f);
    [magh, ~] = read_mag(data);
    s1 = extractAfter(f,'=');
    s2 = extractBefore(s1,')');
    avgH(allfreqs == s2) = mean(magh, "all");
    avgZw(allfreqs == s2) = mean(mage./magh, "all");
end
avgE = 20*log10(avgE);
avgH = 20*log10(avgH);
figure;
plot(freq, avgE, '-s', "LineWidth", 1.5, "Color", "#D95319");
grid on;
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
xlabel("Frequency, f in MHz");
ylabel("Avg. el. field magnitude in dBV/m");
set(gca, "XScale", "log");
xticks([1e-5, 1e-3, 1e-1, 1e0, 1e1, 1e3]);
xticklabels({"10^{-5}", "10^{-3}", "0.1", "1", "10", "10^3"});
figure;
plot(freq, avgH, '-s', "LineWidth", 1.5);
grid on;
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
xlabel("Frequency, f in MHz");
ylabel("Avg. mag. field magnitude in dBA/m");
set(gca, 'XScale', 'Log');
set(gca, "XScale", "log");
xticks([1e-5, 1e-3, 1e-1, 1e0, 1e1, 1e3]);
xticklabels({"10^{-5}", "10^{-3}", "0.1", "1", "10", "10^3"});
figure;
plot(freq, avgZw, '-s', "LineWidth", 1.5, "Color", "#EDB120");
grid on;
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
set(gca, 'XScale', 'Log');
xlabel("Frequency, f in MHz");
ylabel("Avg. wave impedance, Z_W in \Omega");
set(gca, "XScale", "log");
xticks([1e-5, 1e-3, 1e-1, 1e0, 1e1, 1e3]);
xticklabels({"10^{-5}", "10^{-3}", "0.1", "1", "10", "10^3"});
yticks([0, 100, 200, 300, 377, 500]);
ylim([0 500]);
%% Resistor current vs freq
resistor_current = readmatrix("data/fields/resistor/resistor_current.txt");
freq = resistor_current(:,1);
magI = resistor_current(:,2);
figure;
plot(freq, magI, 'LineWidth', 1.5, "Color", "#D95319");
grid on;
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
ylim([0 0.3]);
xlim([0 3e3]);
xlabel('Frequency, f in MHz');
ylabel('Magnitude of current in A');
xticks(0:500:3000);
%% Resistor voltage vs freq
resistor_voltage = readmatrix("data/fields/resistor/resistor_voltage.txt");
freq = resistor_voltage(:,1);
magV = resistor_voltage(:,2);
figure;
plot(freq, magV, 'LineWidth', 1.5);
grid on;
set(gca,'FontSize', 18);
set(gca,'FontName', "Times New Roman");
ylim([0 15]);
xlim([0 3e3]);
xlabel('Frequency, f in MHz');
ylabel("Magnitude of voltage in V");
xticks(0:500:3000);

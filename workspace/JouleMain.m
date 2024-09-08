%% Joule Main

clear;clc;close all;
cd('D:/Works/Analysis/Joule');
addpath(genpath('matnwb'));

generateCore();

disp("Toolbox setup done.");

%%

nwbPath = "data/";
nwbFiles = {dir(nwbPath).name};
nwbFiles = nwbFiles(endsWith(nwbFiles, ".nwb"));
nwbFile = nwbPath + nwbFiles{17};

nwb = nwbRead(nwbFile);
disp("Loaded" + nwbFile);

times = nwb.intervals.get("UNDECODABLE_PRESUMED_CODES").start_time.data(:);
codes = nwb.intervals.get("UNDECODABLE_PRESUMED_CODES").vectordata.get("codes").data(:);
signl = nwb.acquisition.get("probe_1_lfp").electricalseries.get("probe_1_lfp_data").data;
signals = zeros(128, 10000);

for i = 1:length(codes)

    if codes(i) == 4749

        signals = signals + signl(:, times(i)*1000+1:times(i)*1000+10000);

    end

end

%%

a = jSpectrogram(signals);
a(42, :) = (a(43, :) + a(41, :))/2;
a = jMeanFilt2(a, 3, 20);

for i = 1:size(a, 2)

    a(:, i) = a(:, i) / max(a(:, i));

end

imagesc(linspace(0, 500, size(a, 2)), 1:128,log(a));
title("Log scaled power spectrum, scaled to maximum in each frequency across channels");
ylabel("Channel");
xlabel("Frequency");
clim([-1 0]);

%%

cd('D:/Works/Analysis/Joule');
addpath(genpath('matnwb'));

generateCore();
disp("Toolbox setup done.");
nwbPath = "data/";

y3 = jPassiveGLOSignals(nwbPath, 17, "muae", 5, 1, 0, 0, 1, 7);
% y2 = jPassiveGLOSignals(nwbPath, 18, "muae", 5, 1, 1, 1, 1, 7);

% y2 = jPassiveGLOSignals(nwbPath, 3, "muae", 5, 1, 0, 0, 1, 7);
% y3 = jPassiveGLOSignals(nwbPath, 3, "muae", 5, 1, 1, 1, 2, 2);

%%

y = y2;
% y = y2;

%%

save("y.mat", 'y2', 'y1', '-v7.3');
disp("Done");

%%

y = load("y2.mat", 'y2');
disp("Done");

%%

t = 1:6000;
ch1 = 1:128;%3:127;
ch2 = 1:128; % -2
figure("Position", [0 0 1500 1000]);

subplot(3, 1, 1);
img1 = squeeze(mean(y.B2G(:, ch1, t), 1));
imagesc(img1);
hold("on");
colorbar;%clim([-5 5]);
xline(100+5, "LineWidth", .5);
xline(166+5, "LineWidth", .5);
xline(182+5, "LineWidth", .5);
yline(70, "LineWidth", 1.);
title("B2G (Global in GLOEXP)");

subplot(3, 1, 2);
img1 = squeeze(mean(y.B4G(:, ch2, t), 1));clc

imagesc(img1);
hold("on");
colorbar;%clim([-5 5]);
xline(100+5, "LineWidth", 2.5);
xline(166+5, "LineWidth", 2.5);
xline(182+5, "LineWidth", 2.5);
yline(70, "LineWidth", 1.);
title("B4G (Global in sequence control)");

subplot(3, 1, 3);
img2 = squeeze(mean(y.B2G(:, ch1, t), 1)) - squeeze(mean(y.B4G(:, ch2, t), 1));
imagesc(img2);
hold("on");
colorbar;%clim([-2 2]);
xline(100+5, "LineWidth", 2.5);
xline(166+5, "LineWidth", 2.5);
xline(182+5, "LineWidth", 2.5);
yline(70, "LineWidth", 1.);
title("B2G-B4G");

sgtitle(y.Area);

%%
img3 = img2(:, 5035:6000) - img2(:, 4000:4965);
figure("Position", [0 0 1600 900]);
imagesc(img3);
clim([-1.5 1.5]);
title("(B2G-B4G)_{4th Stim} - (B2G-B4G)_{3rd Stim}");
xlabel(y.Area);
% subplot(3, 1, 2);
% phtd = squeeze(mean(y.photodiode_l(:, t), 1));
% plot(t, phtd);

%%

nwbPath = "data/";
nwbFiles = {dir(nwbPath).name};
nwbFiles = nwbFiles(endsWith(nwbFiles, ".nwb"));
nwbFile = nwbPath + nwbFiles{18};

nwb = nwbRead(nwbFile);
% nwb = nwbRead("sub-C31_ses-230303.nwb");

disp("Loaded " + nwbFile);
% GLOPassiveTrialCounts(nwb);

%% RF mapping

close all;
jReceptiveFieldMap(nwb, 2);

%%

tic;
lfps = nwb.acquisition.get("probe_0_lfp").electricalseries.get("probe_0_lfp_data").data(:, :);
muas = nwb.acquisition.get("probe_0_muae").electricalseries.get("probe_0_muae_data").data(:, :);
toc;

%%

close all;
cd('D:/Works/Analysis/Joule');
addpath(genpath('matnwb'));

generateCore();
disp("Toolbox setup done.");
nwbPath = "data/";
y = PassiveGLOSignals(nwbPath, 3, "muae", 5, 2, 0, 1, 1, 7);

%%
figure("Position", [0 0 1500 1000]);
ch = 15;
plot(squeeze(mean(y.B2G(:, ch, 1001:6000), 1)), "DisplayName", "GO");
hold("on");
plot(squeeze(mean(y.B2L(:, ch, 1001:6000), 1)), "DisplayName", "LO");
legend;xlabel("Time(ms)");ylabel("Scaled avg. MUA");
title("Mean MUA single channel [ch." + num2str(ch) + ",~L5/6-Deep-V3d] full-trials");
grid("on");
%%
figure("Position", [0 0 1500 1000]);

chs = 100:110;
subplot(2, 1, 1);
x1 = 501:5000;
y1 = squeeze(mean(mean(y.x1(:, chs, 1501:6000), 1), 2));
e1 = 2*squeeze(mean(std(y.x1(:, chs, 1501:6000), 1), 2))/sqrt(size(y.x1, 1));
errorbar(x1, y1, e1, "DisplayName", "GO+-2SE", "CapSize", 0);

hold("on");

x2 = 501:5000;
y2 = squeeze(mean(mean(y.x2(:, chs, 1501:6000), 1), 2));
e2 = 2*squeeze(mean(std(y.x2(:, chs, 1501:6000), 1), 2))/sqrt(size(y.x2, 1));
errorbar(x2, y2, e2, "DisplayName", "LO+-2SE", "CapSize", 0);

legend;grid;
title("Mean MUA single superficial V3d channels " + num2str(chs(1)) + "-" + num2str(chs(end)) + ", Trial[1-4], errorbar 2*SE");
xlim([500 5000]);

chs = 110:120;
subplot(2, 1, 2);
x1 = 501:5000;
y1 = squeeze(mean(mean(y.x1(:, chs, 1501:6000), 1), 2));
e1 = 2*squeeze(mean(std(y.x1(:, chs, 1501:6000), 1), 2))/sqrt(size(y.x1, 1));
errorbar(x1, y1, e1, "DisplayName", "GO+-2SE", "CapSize", 0);

hold("on");

x2 = 501:5000;
y2 = squeeze(mean(mean(y.x2(:, chs, 1501:6000), 1), 2));
e2 = 2*squeeze(mean(std(y.x2(:, chs, 1501:6000), 1), 2))/sqrt(size(y.x2, 1));
errorbar(x2, y2, e2, "DisplayName", "LO+-2SE", "CapSize", 0);

legend;grid;
title("Mean MUA single deep V3d channels " + num2str(chs(1)) + "-" + num2str(chs(end)) + ", Trial[1-4], errorbar 2*SE");
xlim([500 5000]);
%%
figure("Position", [0 0 1400 1400]);

chs = 5:25;

subplot(3, 1, 1);

x1 = 3001:5000;
y1 = squeeze(mean(mean(y.x1(:, chs, 4001:6000), 1), 2));
% y1 = y1 / max(y1);
e1 = 2*squeeze(mean(std(y.x1(:, chs, 4001:6000), 1), 2))/sqrt(size(y.x1, 1));
errorbar(x1, y1, e1, "DisplayName", "GO+-2SE", "CapSize", 0);

xlabel("Time(ms)");ylabel("Scaled diff.");
hold("on");

x2 = 3001:5000;
y2 = squeeze(mean(mean(y.x2(:, chs, 4001:6000), 1), 2));
% y2 = y2 / max(y1);
e2 = 2*squeeze(mean(std(y.x2(:, chs, 4001:6000), 1), 2))/sqrt(size(y.x2, 1));
errorbar(x2, y2, e2, "DisplayName", "LO+-2SE", "CapSize", 0);

plot(ones(1, 10)*4000, linspace(-.2, .6, 10), "LineWidth", 2, "LineStyle", "--", "DisplayName", "Onset of 4th stim.");

xlabel("Time(ms)");ylabel("Scaled diff.");
legend;grid;
title("V1 mean MUA GLO");
% title("Mean LFP single channel " + num2str(chs(1)) + "-" + num2str(chs(end)) + ", Trial[3-4], errorbar 2*SE");
xlim([3000 5000]);
ylim([-.5 2.5]);

chs = 25:50;

subplot(3, 1, 2);

x1 = 3001:5000;
y1 = squeeze(mean(mean(y.x1(:, chs, 4001:6000), 1), 2));
% y1 = y1 / max(y1);
e1 = 2*squeeze(mean(std(y.x1(:, chs, 4001:6000), 1), 2))/sqrt(size(y.x1, 1));
errorbar(x1, y1, e1, "DisplayName", "GO+-2SE", "CapSize", 0);

xlabel("Time(ms)");ylabel("Scaled diff.");
hold("on");

x2 = 3001:5000;
y2 = squeeze(mean(mean(y.x2(:, chs, 4001:6000), 1), 2));
% y2 = y2 / max(y2);
e2 = 2*squeeze(mean(std(y.x2(:, chs, 4001:6000), 1), 2))/sqrt(size(y.x2, 1));
errorbar(x2, y2, e2, "DisplayName", "LO+-2SE", "CapSize", 0);

plot(ones(1, 10)*4000, linspace(-.5, 2.5, 10), "LineWidth", 2, "LineStyle", "--", "DisplayName", "Onset of 4th stim.");

xlabel("Time(ms)");ylabel("Scaled diff.");
legend;grid;
title("V2 mean MUA GLO");
% title("Mean LFP single channel " + num2str(chs(1)) + "-" + num2str(chs(end)) + ", Trial[3-4], errorbar 2*SE");
xlim([3000 5000]);


chs = 110:128;

subplot(3, 1, 3);

x1 = 3001:5000;
y1 = squeeze(mean(mean(y.x1(:, chs, 4001:6000), 1), 2));
% y1 = y1 / max(y1);
e1 = 2*squeeze(mean(std(y.x1(:, chs, 4001:6000), 1), 2))/sqrt(size(y.x1, 1));
errorbar(x1, y1, e1, "DisplayName", "GO+-2SE", "CapSize", 0);

xlabel("Time(ms)");ylabel("Scaled diff.");
hold("on");

x2 = 3001:5000;
y2 = squeeze(mean(mean(y.x2(:, chs, 4001:6000), 1), 2));
% y2 = y2 / max(y2);
e2 = 2*squeeze(mean(std(y.x2(:, chs, 4001:6000), 1), 2))/sqrt(size(y.x2, 1));
errorbar(x2, y2, e2, "DisplayName", "LO+-2SE", "CapSize", 0);

plot(ones(1, 10)*4000, linspace(-.5, 2.5, 10), "LineWidth", 2, "LineStyle", "--", "DisplayName", "Onset of 4th stim.");

xlabel("Time(ms)");ylabel("Scaled diff.");
legend;grid;
title("V3d mean MUA GLO");
% title("Mean LFP single channel " + num2str(chs(1)) + "-" + num2str(chs(end)) + ", Trial[3-4], errorbar 2*SE");
xlim([3000 5000]);
%% Setup

clear;clc;close all;

cd('D:\Electrophysiology\Matdelane\');
addpath(genpath('matnwb'));
addpath(genpath('matdelane'));
addpath(genpath('flipv2'));
generateCore();

addpath('fieldtrip');ft_defaults;
nwbPath = "data\";
nwbFiles = {dir(nwbPath).name};
nwbFiles = nwbFiles(endsWith(nwbFiles, ".nwb"));
disp("Toolbox setup done.");

%% Ses: sub-C31o_ses-230816 (Probes A,B,C)

%% Probe A (PFC, laminar)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{2};
nwb = nwbRead(nwbFile);
disp(length(nwb.general_extracellular_ephys.keys()));

%% E.1: Load LFP probeA PFC

[c, x] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4500, 0);
% [cm, xm] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4500, 0, "muae");

disp(c{1}.session);

%% E1.1: MUAe plots

imxm = squeeze(mean(xm{12}, 1));
imxm = squeeze(mean(imxm, 1));
imxm = (imxm - mean(imxm)) / std(imxm);
plot(linspace(-500, 4500, 5000), imxm, "DisplayName", "PFC");hold("on");
imxm = squeeze(mean(xm2{12}, 1));
imxm = squeeze(mean(imxm, 1));
imxm = (imxm - mean(imxm)) / std(imxm);
plot(linspace(-500, 4500, 5000), imxm, "DisplayName", "V4/MT");
xlim([-500 4500]);
xline(0, HandleVisibility="off");
xline(1031, HandleVisibility="off");
xline(2062, HandleVisibility="off");
xline(3093, HandleVisibility="off");
title("MUAenv/Zsc/RRRX");
xlabel("Time (ms)");
ylabel("Z-score");
legend;

%% E.2: Channel and layer identification

channel_in_layer = struct();
channel_in_layer.deep = 21:80;
channel_in_layer.mid = 81:85;
channel_in_layer.sup = [86:112, 114:2:128];
goodch = [channel_in_layer.deep, channel_in_layer.mid, channel_in_layer.sup];

jLFPprobeINFO(x{1}(:, goodch, :));

%% E.3: Evaluate vFLIP

jVFLIP(x{1}(:, goodch, :));

%% E.4: TFR calculations all trials; PFC

freqlims = [0 200];
leakage = 0.85;
freqres = 5.0;
overlap = 95;

y = squeeze(mean(mean(x{1}, 1), 2));
[p1, f1, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", freqlims, "OverlapPercent", overlap, "FrequencyResolution", freqres, "Leakage", leakage);
pgx = cell(12, 3);

for ik = [1, 5, 9]
    
    xG1d = squeeze(mean(x{ik}(:, channel_in_layer.deep, :), 2));
    xG1m = squeeze(mean(x{ik}(:, channel_in_layer.mid, :), 2));
    xG1s = squeeze(mean(x{ik}(:, channel_in_layer.sup, :), 2));
    
    TR = size(xG1d, 1);

    pgxcd = zeros([TR, size(p1)]);
    pgxcm = zeros([TR, size(p1)]);
    pgxcs = zeros([TR, size(p1)]);
    
    parfor jk = 1:TR % trials
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1d(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
        pgxcd(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
        pgxcm(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
        pgxcs(jk, :, :) = p1temp;
    
        if mod(jk, 20) == 0

            fprintf(num2str(jk));

        end
    
    end

    pgx{ik, 1} = pgxcd;
    pgx{ik, 2} = pgxcm;
    pgx{ik, 3} = pgxcs;

    disp(" >cond : " + num2str(ik));

end

%% E.5: Band and Time specifications

fmap = f1;
tmap = (t1 - 0.500 - t1(1) + 0.031)*1000;% TFR Kaiser's time window offset shift = t(1)*2

tbands = cell(1, 5);
tbandlabels = ["Fix", "S1d1", "S2d2", "S3d3", "S4d4"];
tbands{1} = 1:find(tmap > 0, 1);
tbands{2} = find(tmap > 0, 1):find(tmap > 1000, 1);
tbands{3} = find(tmap > 1031, 1):find(tmap > 2031, 1);
tbands{4} = find(tmap > 2062, 1):find(tmap > 3062, 1);
tbands{5} = find(tmap > 3093, 1):find(tmap > 4093, 1);

fbands = cell(1, 5);
fbandlabels = ["Theta(2-7Hz)", "Alpha(8-12Hz)", "Beta(14-30Hz)", "GammaL(32-80Hz)", "GammaH(80+Hz)"];
fbands{1} = find(fmap > 2, 1):find(fmap > 7, 1);
fbands{2} = find(fmap > 8, 1):find(fmap > 12, 1);
fbands{3} = find(fmap > 14, 1):find(fmap > 30, 1);
fbands{4} = find(fmap > 32, 1):find(fmap > 80, 1);
fbands{5} = find(fmap > 80, 1):find(fmap >= max(fmap), 1);

%% E.6: Visualize TFR

condinflabel = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA",...
    "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

layeridlabel = ["deep", "mid", "sup"];
areainf = "PFC/";

tcond1 = 9;
layerid = 1;

figure;
subplot(2, 1, 1);
tfr1 = squeeze(mean(pgx{tcond1, layerid}, 1));

for ik = 1:size(tfr1, 1)

    tfr1(ik, :) = tfr1(ik, :) / mean(tfr1(ik, tbands{1}));

end

imagesc(10*log10(tfr1), "XData", tmap, "YData", fmap);
% ylim([0 20])
xlim([tmap(1) tmap(end)]);
set(gca, "YDir", "normal");
hold("on");
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time (ms)");
ylabel("Frequency (Hz)");
title("Power change from baseline (dB)");
colorbar;

for fband = 1:5

    yline(fmap(fbands{fband}(1)), "Color", [1 0 0]);

end
xlabel("Time (ms)");

subplot(2, 1, 2);

for fband = 1:5

    tfr1 = squeeze(mean(pgx{tcond1, layerid}(:, fbands{fband}, :), 1));
    tfr1 = squeeze(mean(tfr1, 1));
    tfr1 = tfr1 / mean(tfr1(tbands{1}));
    plot(tmap, 10*log10(tfr1), "DisplayName", fbandlabels(fband), "LineWidth", 2);
    % ylim([0 20])
    xlim([tmap(1) tmap(end)]);
    hold("on");
    xline(0, HandleVisibility="off");
    xline(1031, HandleVisibility="off");
    xline(2062, HandleVisibility="off");
    xline(3093, HandleVisibility="off");
    xlabel("Times (ms)");
    ylabel("Power change (dB)");
    title("Power change from fixation baseline");

end

colorbar;
sgtitle(areainf + condinflabel(tcond1) + "/" + layeridlabel(layerid));
legend;

%% E.7: Band PEV

layerid = 1;
condinf = [1, 5];

layerinf = layeridlabel(layerid) + " layer";
expvars = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, :);
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, :);
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    
    data = zeros(N1+N2, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    % data = jSmooth(data, 50);
    
    groupIDs = [ones(1, N1), ones(1, N2)*2];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars{fband} = squeeze(expv)*100;

end

%% E.8: PEV plot

figure;

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars{fband}) / sqrt(size(expvars{fband}, 1));
    plot(tmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(tmap, stx, "Color", cl);
    plot(tmap, sty, "Color", cl);
    patch([tmap', tmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(0);
    xline(1031);
    xline(2062);
    xline(3093);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " " + condinflabel(condinf(1)) + "-vs-" + condinflabel(condinf(2)) + " PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf);


%% Probe B (V4/MT, V4 laminar, MT partial coverage)

%% E.1: Load LFP probeB  V4-MT

[c, x] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4500, 1);
disp(c{1}.session);

%% E.2: Channel and layer identification


channel_in_layer = struct();
channel_in_layer.deep = 1:27;
channel_in_layer.mid = 29:2:35;
channel_in_layer.sup = [36:44, 46:64];
goodch = [channel_in_layer.deep, channel_in_layer.mid, channel_in_layer.sup];

channel_in_layer2 = struct();
channel_in_layer2.deep = 81:128;
channel_in_layer2.mid = 76:80;
channel_in_layer2.sup = 66:75;
goodch2 = [channel_in_layer2.sup, channel_in_layer2.mid, channel_in_layer2.deep];

jLFPprobeINFO(x{1}(:, goodch, :));
jLFPprobeINFO(x{1}(:, goodch2, :));

%% E.3: Evaluate vFLIP

jVFLIP(x{1}(:, goodch, :));
jVFLIP(x{1}(:, goodch2, :));

%% E.4: TFR calculations all trials; V4

freqlims = [0 200];
leakage = 0.85;
freqres = 5.0;
overlap = 95;

y = squeeze(mean(mean(x{1}, 1), 2));
[p1, f1, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", freqlims, "OverlapPercent", overlap, "FrequencyResolution", freqres, "Leakage", leakage);
pgx = cell(12, 3);

for ik = [1, 5, 9]
    
    xG1d = squeeze(mean(x{ik}(:, channel_in_layer.deep, :), 2));
    xG1m = squeeze(mean(x{ik}(:, channel_in_layer.mid, :), 2));
    xG1s = squeeze(mean(x{ik}(:, channel_in_layer.sup, :), 2));
    
    TR = size(xG1d, 1);

    pgxcd = zeros([TR, size(p1)]);
    pgxcm = zeros([TR, size(p1)]);
    pgxcs = zeros([TR, size(p1)]);
    
    parfor jk = 1:TR % trials
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1d(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
        pgxcd(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
        pgxcm(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
        pgxcs(jk, :, :) = p1temp;
    
        if mod(jk, 20) == 0

            fprintf(num2str(jk));

        end
    
    end

    pgx{ik, 1} = pgxcd;
    pgx{ik, 2} = pgxcm;
    pgx{ik, 3} = pgxcs;

    disp(" >cond : " + num2str(ik));

end

%% E.5: Band and Time specifications

fmap = f1;
tmap = (t1 - 0.500 - t1(1) + 0.031)*1000;% TFR Kaiser's time window offset shift = t(1)*2

tbands = cell(1, 5);
tbandlabels = ["Fix", "S1d1", "S2d2", "S3d3", "S4d4"];
tbands{1} = 1:find(tmap > 0, 1);
tbands{2} = find(tmap > 0, 1):find(tmap > 1000, 1);
tbands{3} = find(tmap > 1031, 1):find(tmap > 2031, 1);
tbands{4} = find(tmap > 2062, 1):find(tmap > 3062, 1);
tbands{5} = find(tmap > 3093, 1):find(tmap > 4093, 1);

fbands = cell(1, 5);
fbandlabels = ["Theta(2-7Hz)", "Alpha(8-12Hz)", "Beta(14-30Hz)", "GammaL(32-80Hz)", "GammaH(80+Hz)"];
fbands{1} = find(fmap > 2, 1):find(fmap > 7, 1);
fbands{2} = find(fmap > 8, 1):find(fmap > 12, 1);
fbands{3} = find(fmap > 14, 1):find(fmap > 30, 1);
fbands{4} = find(fmap > 32, 1):find(fmap > 80, 1);
fbands{5} = find(fmap > 80, 1):find(fmap >= max(fmap), 1);

%% E.6: Visualize TFR

condinflabel = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA",...
    "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

layeridlabel = ["deep", "mid", "sup"];
areainf = "V4/";

tcond1 = 12;
layerid = 1;

figure;
subplot(2, 1, 1);
tfr1 = squeeze(mean(pgx{tcond1, layerid}, 1));

for ik = 1:size(tfr1, 1)

    tfr1(ik, :) = tfr1(ik, :) / mean(tfr1(ik, tbands{1}));

end

imagesc(10*log10(tfr1), "XData", tmap, "YData", fmap);
% ylim([0 20])
xlim([tmap(1) tmap(end)]);
set(gca, "YDir", "normal");
hold("on");
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time (ms)");
ylabel("Frequency (Hz)");
title("Power change from baseline (dB)");
colorbar;

for fband = 1:5

    yline(fmap(fbands{fband}(1)), "Color", [1 0 0]);

end
xlabel("Time (ms)");

subplot(2, 1, 2);

for fband = 1:5

    tfr1 = squeeze(mean(pgx{tcond1, layerid}(:, fbands{fband}, :), 1));
    tfr1 = squeeze(mean(tfr1, 1));
    tfr1 = tfr1 / mean(tfr1(tbands{1}));
    plot(tmap, 10*log10(tfr1), "DisplayName", fbandlabels(fband), "LineWidth", 2);
    % ylim([0 20])
    xlim([tmap(1) tmap(end)]);
    hold("on");
    xline(0, HandleVisibility="off");
    xline(1031, HandleVisibility="off");
    xline(2062, HandleVisibility="off");
    xline(3093, HandleVisibility="off");
    xlabel("Times (ms)");
    ylabel("Power change (dB)");
    title("Power change from fixation baseline");

end

colorbar;
sgtitle(areainf + condinflabel(tcond1) + "/" + layeridlabel(layerid));
legend;

%% E.7: Band PEV

layerid = 1;
condinf = [4, 8];

layerinf = layeridlabel(layerid) + " layer";
expvars = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, :);
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, :);
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    
    data = zeros(N1+N2, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    % data = jSmooth(data, 50);
    
    groupIDs = [ones(1, N1), ones(1, N2)*2];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars{fband} = squeeze(expv)*100;

end

%% E.8: PEV plot

figure;

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars{fband}) / sqrt(size(expvars{fband}, 1));
    plot(tmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(tmap, stx, "Color", cl);
    plot(tmap, sty, "Color", cl);
    patch([tmap', tmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(0);
    xline(1031);
    xline(2062);
    xline(3093);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " " + condinflabel(condinf(1)) + "-vs-" + condinflabel(condinf(2)) + " PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf);

%%
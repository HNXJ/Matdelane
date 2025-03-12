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

%% Ses: sub-V198o_ses-230721 (Probes A,B)

%% Probe A (PFC, laminar)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{16};
nwb = nwbRead(nwbFile);
disp(length(nwb.general_extracellular_ephys.keys()));

layeridlabel = ["deep", "mid", "sup"];
areainf = "PFC/";

condinflabel = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA",...
    "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

%% E.1: Load LFP probeA PFC

[c, x] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4000, 0);
disp(c{1}.session);

%% E1.1: MUAe plots

[cm, xm] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4000, 0, "muae");
disp(cm{1}.session);

condid = 12;
figure;

imxm = squeeze(mean(xm{condid}, 1));
imxm = squeeze(mean(imxm, 1));
imxm = (imxm - mean(imxm)) / std(imxm);
plot(linspace(-500, 4000, 5000), imxm, "DisplayName", areainf);

hold("on");
xline(0, HandleVisibility="off");
xline(1031, HandleVisibility="off");
xline(2062, HandleVisibility="off");
xline(3093, HandleVisibility="off");

title("MUAenv/Zsc/" + condinflabel(condid));
xlabel("Time (ms)");
ylabel("Z-score");
xlim([-500 4000]);

legend;

%% E.2: Channel and layer identification

channel_in_layer = struct();
channel_in_layer.deep = 21:81;
channel_in_layer.mid = 82:86;
channel_in_layer.sup = [87:112, 114:2:128];
goodch = [channel_in_layer.deep, channel_in_layer.mid, channel_in_layer.sup];

jLFPprobeINFO(x{1}(:, goodch, :));

%% E.3: Evaluate vFLIP

jVFLIP(x{1}(:, goodch, :));

%% E.4: TFR calculations all trials; PFC

freqlims = [0 200];
leakage = 0.85;
freqres = 5.0;
overlap = 95;

channel_in_layer_selected = channel_in_layer;
y = squeeze(mean(mean(x{1}, 1), 2));
[p1, f1, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", freqlims, "OverlapPercent", overlap, "FrequencyResolution", freqres, "Leakage", leakage);
pgx = cell(12, 3);

for ik = 1:12
    
    xG1d = squeeze(mean(x{ik}(:, channel_in_layer_selected.deep, :), 2));
    xG1m = squeeze(mean(x{ik}(:, channel_in_layer_selected.mid, :), 2));
    xG1s = squeeze(mean(x{ik}(:, channel_in_layer_selected.sup, :), 2));
    
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
tmap = (t1 - 0.500)*1000;% TFR Kaiser's time window offset shift = t(1)*2

tbands = cell(1, 5);
tbandlabels = ["Fix", "S1d1", "S2d2", "S3d3", "S4d4"];

tbands{1} = 1:find(tmap > -50, 1);
tbands{2} = find(tmap > 0, 1):find(tmap > 1000, 1);
tbands{3} = find(tmap > 1031, 1):find(tmap > 2031, 1);
tbands{4} = find(tmap > 2062, 1):find(tmap > 3062, 1);
tbands{5} = find(tmap > 3093, 1):find(tmap > 4093, 1);

nt_temp = length(tbands{3});

for ik = 2:5

    tbands{ik} = tbands{ik}(1:nt_temp);

end

fbands = cell(1, 5);
fbandlabels = ["Theta(2-7Hz)", "Alpha(8-12Hz)", "Beta(14-30Hz)", "GammaL(32-80Hz)", "GammaH(80+Hz)"];

fbands{1} = find(fmap > 2, 1):find(fmap > 7, 1);
fbands{2} = find(fmap > 8, 1):find(fmap > 12, 1);
fbands{3} = find(fmap > 14, 1):find(fmap > 30, 1);
fbands{4} = find(fmap > 32, 1):find(fmap > 80, 1);
fbands{5} = find(fmap > 80, 1):find(fmap >= max(fmap), 1);

%% E.6: Visualize TFR

tcond1 = 1;
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

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 1;
condinf = [2, 3, 4];

layerinf2 = layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars2{fband} = squeeze(expv)*100;

end

%% E.10: Omission PEV (Positional, BX2/BX3/BX4)

layerid = 1;
condinf = [6, 7, 8];

layerinf3 = layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars3{fband} = squeeze(expv)*100;

end

%% E.11: Omission PEV (Positional, RX2/RX3/RX4)

layerid = 1;
condinf = [10, 11, 12];

layerinf4 = layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars4{fband} = squeeze(expv)*100;

end

%% E.12: PEV plot AX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars2{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars2{fband}) / sqrt(size(expvars2{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% E.13: PEV plot BX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars3{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars3{fband}) / sqrt(size(expvars3{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Bx/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% E.14: PEV plot RX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars4{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars4{fband}) / sqrt(size(expvars4{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Rx/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% Probe B (V4/MT, It seems to be majorly MT due to strong gamma response)

%% E.1: Load LFP probeB  V4-MT

[c, x] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4000, 1);
disp(c{1}.session);
areainf = "V4/MT";

%% E1.1: MUAe plots

[cm2, xm2] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4000, 1, "muae");
disp(cm2{1}.session);

condid = 12;
figure;

imxm = squeeze(mean(xm2{condid}, 1));
imxm = squeeze(mean(imxm, 1));
imxm = (imxm - mean(imxm)) / std(imxm);
plot(linspace(-500, 4000, 5000), imxm, "DisplayName", areainf);

hold("on");
xline(0, HandleVisibility="off");
xline(1031, HandleVisibility="off");
xline(2062, HandleVisibility="off");
xline(3093, HandleVisibility="off");

title("MUAenv/Zsc/" + condinflabel(condid));
xlabel("Time (ms)");
ylabel("Z-score");
xlim([-500 4000]);

legend;

%% E.2: Channel and layer identification

channel_in_layer = struct();
channel_in_layer.deep = 1:27;
channel_in_layer.mid = 29:2:35;
channel_in_layer.sup = [36:44, 46:60];
goodch = [channel_in_layer.deep, channel_in_layer.mid, channel_in_layer.sup];

channel_in_layer2 = struct();
channel_in_layer2.deep = 81:128;
channel_in_layer2.mid = 76:80;
channel_in_layer2.sup = [61:64, 66:75];
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

channel_in_layer_selected = channel_in_layer;
y = squeeze(mean(mean(x{1}, 1), 2));
[p1, f1, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", freqlims, "OverlapPercent", overlap, "FrequencyResolution", freqres, "Leakage", leakage);
pgx = cell(12, 3);

for ik = 1:12
    
    xG1d = squeeze(mean(x{ik}(:, channel_in_layer_selected.deep, :), 2));
    xG1m = squeeze(mean(x{ik}(:, channel_in_layer_selected.mid, :), 2));
    xG1s = squeeze(mean(x{ik}(:, channel_in_layer_selected.sup, :), 2));
    
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
tmap = (t1 - 0.500)*1000;% TFR Kaiser's time window offset shift = t(1)*2

tbands = cell(1, 5);
tbandlabels = ["Fix", "S1d1", "S2d2", "S3d3", "S4d4"];

tbands{1} = 1:find(tmap > -50, 1);
tbands{2} = find(tmap > 0, 1):find(tmap > 1000, 1);
tbands{3} = find(tmap > 1031, 1):find(tmap > 2031, 1);
tbands{4} = find(tmap > 2062, 1):find(tmap > 3062, 1);
tbands{5} = find(tmap > 3093, 1):find(tmap > 4093, 1);

nt_temp = length(tbands{3});

for ik = 2:5

    tbands{ik} = tbands{ik}(1:nt_temp);

end

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

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 1;
condinf = [2, 3, 4];

layerinf2 = layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars2{fband} = squeeze(expv)*100;

end

%% E.10: Omission PEV (Positional, BX2/BX3/BX4)

layerid = 1;
condinf = [6, 7, 8];

layerinf3 = layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars3{fband} = squeeze(expv)*100;

end

%% E.11: Omission PEV (Positional, RX2/RX3/RX4)

layerid = 1;
condinf = [10, 11, 12];

layerinf4 = layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars4{fband} = squeeze(expv)*100;

end

%% E.12: PEV plot AX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars2{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars2{fband}) / sqrt(size(expvars2{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% E.13: PEV plot BX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars3{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars3{fband}) / sqrt(size(expvars3{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Bx/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% E.14: PEV plot RX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars4{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars4{fband}) / sqrt(size(expvars4{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Rx/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% Probe C (V1/V2, It seems that it consists only of one area, to be identified, highly foveal)

%% E.1: Load LFP probeC  V1

[c, x] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4000, 2);
disp(c{1}.session);
areainf = "V1";

%% E1.1: MUAe plots

[cm2, xm2] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4000, 2, "muae");
disp(cm2{1}.session);

condid = 12;
figure;

imxm = squeeze(mean(xm2{condid}, 1));
imxm = squeeze(mean(imxm, 1));
imxm = (imxm - mean(imxm)) / std(imxm);
plot(linspace(-500, 4000, 5000), imxm, "DisplayName", areainf);

hold("on");
xline(0, HandleVisibility="off");
xline(1031, HandleVisibility="off");
xline(2062, HandleVisibility="off");
xline(3093, HandleVisibility="off");

title("MUAenv/Zsc/" + condinflabel(condid));
xlabel("Time (ms)");
ylabel("Z-score");
xlim([-500 4000]);

legend;

%% E.2: Channel and layer identification

channel_in_layer = struct(); % V1
channel_in_layer.deep = [46:107, 109:115];
channel_in_layer.mid = 41:45;
channel_in_layer.sup = 1:40;
goodch = [channel_in_layer.sup, channel_in_layer.mid, channel_in_layer.deep];

% channel_in_layer2 = struct(); % V2
% channel_in_layer2.sup = 81:100;
% channel_in_layer2.mid = 101:105;
% channel_in_layer2.deep = [106:107, 109:128];
% goodch2 = [channel_in_layer2.sup, channel_in_layer2.mid, channel_in_layer2.deep];

jLFPprobeINFO(x{1}(:, goodch, :));
% jLFPprobeINFO(x{1}(:, goodch2, :));

%% E.3: Evaluate vFLIP

jVFLIP(x{1}(:, goodch, :));
% jVFLIP(x{1}(:, goodch2, :));

%% E.4: TFR calculations all trials; V1

freqlims = [0 200];
leakage = 0.85;
freqres = 5.0;
overlap = 95;

channel_in_layer_selected = channel_in_layer;
y = squeeze(mean(mean(x{1}, 1), 2));
[p1, f1, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", freqlims, "OverlapPercent", overlap, "FrequencyResolution", freqres, "Leakage", leakage);
pgx = cell(12, 3);

for ik = 1:12
    
    xG1d = squeeze(mean(x{ik}(:, channel_in_layer_selected.deep, :), 2));
    xG1m = squeeze(mean(x{ik}(:, channel_in_layer_selected.mid, :), 2));
    xG1s = squeeze(mean(x{ik}(:, channel_in_layer_selected.sup, :), 2));
    
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
tmap = (t1 - 0.500)*1000;% TFR Kaiser's time window offset shift = t(1)*2

tbands = cell(1, 5);
tbandlabels = ["Fix", "S1d1", "S2d2", "S3d3", "S4d4"];

tbands{1} = 1:find(tmap > -50, 1);
tbands{2} = find(tmap > 0, 1):find(tmap > 1000, 1);
tbands{3} = find(tmap > 1031, 1):find(tmap > 2031, 1);
tbands{4} = find(tmap > 2062, 1):find(tmap > 3062, 1);
tbands{5} = find(tmap > 3093, 1):find(tmap > 4093, 1);

nt_temp = length(tbands{3});

for ik = 2:5

    tbands{ik} = tbands{ik}(1:nt_temp);

end

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

tcond1 = 8;
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

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 1;
condinf = [2, 3, 4];

layerinf2 = layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars2{fband} = squeeze(expv)*100;

end

%% E.10: Omission PEV (Positional, BX2/BX3/BX4)

layerid = 1;
condinf = [6, 7, 8];

layerinf3 = layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars3{fband} = squeeze(expv)*100;

end

%% E.11: Omission PEV (Positional, RX2/RX3/RX4)

layerid = 1;
condinf = [10, 11, 12];

layerinf4 = layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = pgx{condinf(1), layerid}(:, fbands{fband}, tbands{3}); % S2d2
    x2 = pgx{condinf(2), layerid}(:, fbands{fband}, tbands{4}); % S3d3
    x3 = pgx{condinf(3), layerid}(:, fbands{fband}, tbands{5}); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    expvars4{fband} = squeeze(expv)*100;

end

%% E.12: PEV plot AX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars2{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars2{fband}) / sqrt(size(expvars2{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% E.13: PEV plot BX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars3{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars3{fband}) / sqrt(size(expvars3{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Bx/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%% E.14: PEV plot RX

figure;
ntmap = tmap(tbands{2});

for fband = 1:5

    subplot(5, 1, fband);
    yt = mean(expvars4{fband}, 1);
    yt = yt - min(yt);
    st = std(expvars4{fband}) / sqrt(size(expvars4{fband}, 1));
    plot(ntmap, yt);hold("on");
    stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
    sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
    cl = [0.9 0.7 0.7];
    plot(ntmap, stx, "Color", cl);
    plot(ntmap, sty, "Color", cl);
    patch([ntmap', ntmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
    xline(500);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(fbandlabels(fband));

end

sgtitle("Area:" + areainf + " posOmission/Rx/PEV/TFR/+-2SEM/fRes=" + num2str(freqres) + "Hz/ovlrp=." + num2str(overlap) + " " + layerinf2);

%%
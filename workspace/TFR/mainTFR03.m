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

%% Ses: sub-C31o_ses-230818 (Probes A,B,C)

%% Probe A (PFC, laminar)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{3};
nwb = nwbRead(nwbFile);
disp(length(nwb.general_extracellular_ephys.keys()));

%% E.1: Load LFP probeA PFC

[c, x] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4500, 0);
disp(c{1}.session);

%% E.3: Evaluate channels

deepch_s3_mt = 31:40;
midch_s3_mt = 41:45;
supch_s3_mt = 46:75;

supch_s3_mst = 78:99;
midch_s3_mst = 100:107;
deepch_s3_mst = 109:128;

goodch_mt = [deepch_s3_mt, midch_s3_mt, supch_s3_mt];
goodch_mst = [supch_s3_mst, midch_s3_mst, deepch_s3_mst];

% jLFPprobeINFO(x{1}(:, :, :));
jLFPprobeINFO(x{1}(:, goodch_mst, :));

%% E.4: Evaluate vFLIP

% jVFLIP(x{1}(:, :, :));
jVFLIP(x{1}(:, goodch_mst, :));

%% E.5: TFR calculations

freqlims = [0 200];
timeres = 0.1;
overlap = 95;

y = squeeze(mean(mean(x{1}, 1), 2));
[p1, f1, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", freqlims, "OverlapPercent", overlap, "FrequencyResolution", 20);
pgx = cell(12, 3);

for ik = [1, 5]
    
    xG1d = squeeze(mean(x{ik}(:, deepch_s3_mst, :), 2));
    xG1m = squeeze(mean(x{ik}(:, midch_s3_mst, :), 2));
    xG1s = squeeze(mean(x{ik}(:, supch_s3_mst, :), 2));
    
    TR = size(xG1d, 1);

    pgxcd = zeros([TR, size(p1)]);
    pgxcm = zeros([TR, size(p1)]);
    pgxcs = zeros([TR, size(p1)]);
    
    for jk = 1:TR % trials
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1d(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", 20, "OverlapPercent", overlap);
        pgxcd(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", 20, "OverlapPercent", overlap);
        pgxcm(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", 20, "OverlapPercent", overlap);
        pgxcs(jk, :, :) = p1temp;
    
        if mod(jk, 10) == 0

            fprintf(num2str(jk));

        end
    
    end

    pgx{ik, 1} = pgxcd;
    pgx{ik, 2} = pgxcm;
    pgx{ik, 3} = pgxcs;

    disp("cond : " + num2str(ik));

end

%% E.6: Visualize TFR

figure;
tfr1 = squeeze(mean(pgx{1, 1}, 1));
imagesc(log10(tfr1));

%% E.7: Band specifications

fcnt = size(p1, 1);
tcnt = size(p1, 2);
fmap = linspace(0, 200, fcnt);
tmap = linspace(-500, 4500, tcnt);

thetamap = find(fmap > 2, 1):find(fmap > 7, 1);
alphamap = find(fmap > 8, 1):find(fmap > 12, 1);
betamap = find(fmap > 14, 1):find(fmap > 30, 1);
gammamap1 = find(fmap > 32, 1):find(fmap > 64, 1);
gammamap2 = find(fmap > 64, 1):find(fmap >= 200, 1);

%% E.9: Band PEV

condinflabel = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA",...
    "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

layeridlabel = ["deep", "mid", "sup"];
areainf = "MST : ";
layerid = 3;
condinf = [1, 5];

layerinf = layeridlabel(layerid) + " layer ";

x1 = pgx{condinf(1), layerid}(:, thetamap, :);
x2 = pgx{condinf(2), layerid}(:, thetamap, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars0 = squeeze(expv)*100;

x1 = pgx{condinf(1), layerid}(:, alphamap, :);
x2 = pgx{condinf(2), layerid}(:, alphamap, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars1 = squeeze(expv)*100;

x1 = pgx{condinf(1), layerid}(:, betamap, :);
x2 = pgx{condinf(2), layerid}(:, betamap, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars2 = squeeze(expv)*100;

x1 = pgx{condinf(1), layerid}(:, gammamap2, :);
x2 = pgx{condinf(2), layerid}(:, gammamap2, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars3 = squeeze(expv)*100;

% E.9: PEV plot

figure;

subplot(4, 1, 1);
yt = mean(expvars0, 1);
st = std(expvars0) / sqrt(size(expvars0, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Theta[2-7Hz]");

subplot(4, 1, 2);
yt = mean(expvars1, 1);
st = std(expvars1) / sqrt(size(expvars1, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Alpha[8-12Hz]");

subplot(4, 1, 3);
yt = mean(expvars2, 1);
st = std(expvars2) / sqrt(size(expvars2, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Beta[13-30Hz]");

subplot(4, 1, 4);
yt = mean(expvars3, 1);
st = std(expvars3) / sqrt(size(expvars3, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Gamma[32+Hz]");

sgtitle("Area " + areainf + " " + condinflabel(condinf(1)) + " vs " + condinflabel(condinf(2)) + " PEV in TFR bands +- 2SEM fres=1Hz, ovlrp=.95 " + layerinf);

%% E.2: Load LFP ses 3/17 probeA PFC; Agranular

[c, x] = jOGLOSignals(nwb, "omission_glo_passive", 500, 4500, 0);

%% E.3: Evaluate channels

deepch_s3_pfc = 31:40;
midch_s3_pfc = 41:45;
supch_s3_pfc = 46:75;

goodch_pfc = [deepch_s3_pfc, midch_s3_pfc, supch_s3_pfc];

% jLFPprobeINFO(x{1}(:, :, :));
jLFPprobeINFO(x{1}(:, goodch_pfc, :));

%% E.4: Evaluate vFLIP

jVFLIP(x{1}(:, 22:128, :));
% jVFLIP(x{1}(:, goodch_pfc, :));

%% E.5: TFR calculations

freqlims = [0 200];
timeres = 0.1;
overlap = 95;

y = squeeze(mean(mean(x{1}, 1), 2));
[p1, f1, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", freqlims, "OverlapPercent", overlap, "FrequencyResolution", 20);
pgx = cell(12, 3);

for ik = [1, 5]
    
    xG1d = squeeze(mean(x{ik}(:, deepch_s3_mst, :), 2));
    xG1m = squeeze(mean(x{ik}(:, midch_s3_mst, :), 2));
    xG1s = squeeze(mean(x{ik}(:, supch_s3_mst, :), 2));
    
    TR = size(xG1d, 1);

    pgxcd = zeros([TR, size(p1)]);
    pgxcm = zeros([TR, size(p1)]);
    pgxcs = zeros([TR, size(p1)]);
    
    for jk = 1:TR % trials
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1d(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", 20, "OverlapPercent", overlap);
        pgxcd(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", 20, "OverlapPercent", overlap);
        pgxcm(jk, :, :) = p1temp;
    
        [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", 20, "OverlapPercent", overlap);
        pgxcs(jk, :, :) = p1temp;
    
        if mod(jk, 10) == 0

            fprintf(num2str(jk));

        end
    
    end

    pgx{ik, 1} = pgxcd;
    pgx{ik, 2} = pgxcm;
    pgx{ik, 3} = pgxcs;

    disp("cond : " + num2str(ik));

end

%% E.6: Visualize TFR

figure;
tfr1 = squeeze(mean(pgx{1, 1}, 1));
imagesc(log10(tfr1));

%% E.7: Band specifications

fcnt = size(p1, 1);
tcnt = size(p1, 2);
fmap = linspace(0, 200, fcnt);
tmap = linspace(-500, 4500, tcnt);

thetamap = find(fmap > 2, 1):find(fmap > 7, 1);
alphamap = find(fmap > 8, 1):find(fmap > 12, 1);
betamap = find(fmap > 14, 1):find(fmap > 30, 1);
gammamap1 = find(fmap > 32, 1):find(fmap > 64, 1);
gammamap2 = find(fmap > 64, 1):find(fmap >= 200, 1);

%% E.9: Band PEV

condinflabel = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA",...
    "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

layeridlabel = ["deep", "mid", "sup"];
areainf = "MST : ";
layerid = 3;
condinf = [1, 5];

layerinf = layeridlabel(layerid) + " layer ";

x1 = pgx{condinf(1), layerid}(:, thetamap, :);
x2 = pgx{condinf(2), layerid}(:, thetamap, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars0 = squeeze(expv)*100;

x1 = pgx{condinf(1), layerid}(:, alphamap, :);
x2 = pgx{condinf(2), layerid}(:, alphamap, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars1 = squeeze(expv)*100;

x1 = pgx{condinf(1), layerid}(:, betamap, :);
x2 = pgx{condinf(2), layerid}(:, betamap, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars2 = squeeze(expv)*100;

x1 = pgx{condinf(1), layerid}(:, gammamap2, :);
x2 = pgx{condinf(2), layerid}(:, gammamap2, :);

[N1, nF, nT] = size(x1);
N2 = size(x2, 1);

data = zeros(N1+N2, nF, nT);
data(1:N1, :, :) = x1;
data(N1+1:N1+N2, :, :) = x2;
% data = jSmooth(data, 50);

groupIDs = [ones(1, N1), ones(1, N2)*2];
[expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
expvars3 = squeeze(expv)*100;

% E.9: PEV plot

figure;

subplot(4, 1, 1);
yt = mean(expvars0, 1);
st = std(expvars0) / sqrt(size(expvars0, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Theta[2-7Hz]");

subplot(4, 1, 2);
yt = mean(expvars1, 1);
st = std(expvars1) / sqrt(size(expvars1, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Alpha[8-12Hz]");

subplot(4, 1, 3);
yt = mean(expvars2, 1);
st = std(expvars2) / sqrt(size(expvars2, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Beta[13-30Hz]");

subplot(4, 1, 4);
yt = mean(expvars3, 1);
st = std(expvars3) / sqrt(size(expvars3, 1));
plot(tmap, yt);hold("on");
stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
cl = [0.9 0.7 0.7];
plot(tmap, stx, "Color", cl);
plot(tmap, sty, "Color", cl);
patch([tmap, tmap(end:-1:1)], [stx;sty(end:-1:1)], cl);
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time(ms)");ylabel("PEV(%)");
title("Gamma[32+Hz]");

sgtitle("Area " + areainf + " " + condinflabel(condinf(1)) + " vs " + condinflabel(condinf(2)) + " PEV in TFR bands +- 2SEM fres=1Hz, ovlrp=.95 " + layerinf);

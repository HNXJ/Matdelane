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

%% E.0.1: jNWB object

q1 = jnwb(nwbFile, "PFC/", 500, 4250, 0, 0);

%% E.0.2: MUA plot

q1.jMUAplot(9, [-500 4000]);

%% E.0.3: SUA plot

q1.jSUAplot(9, [100 4000], 100:120);

%% E.0.4: Rastrogram

q1.jRastrogram(1, [1000 3000], 10:120, 1:2);

%% E.1: Channel and layer specs

channel_in_layer = struct();
channel_in_layer.deep = [21:65, 67:81];
channel_in_layer.mid = 82:86;
channel_in_layer.sup = [87:111, 112:2:128];
channel_in_layer.goodch = [channel_in_layer.deep, channel_in_layer.mid, channel_in_layer.sup];

q1.channelinfo{1} = channel_in_layer;

%% E.2: LFP info plot

q1.jLFPprobeINFO(channel_in_layer.goodch, 3);

%% E.3: Evaluate vFLIP

a1 = q1.jVFLIP(channel_in_layer.goodch, 3601:4200, 12);
a2 = q1.jVFLIP(channel_in_layer.goodch, 4101:4700, 12);

imx1 = a1.relpow - a2.relpow;
imx1 = smoothdata2(imx1, "movmedian", "omitmissing", "SmoothingFactor", 0.2);
figure;
imagesc(imx1, "XData", linspace(0, 150, size(imx1, 1)));
clim([-.25 .25]);
colorbar;
xlabel("Frequency");
ylabel("Channel");
title("Omission - Baseline");

%% E.4: TFR calculations all trials

% q1.jCalcTFRs(channel_in_layer);
q1.jCalcTFRs(channel_in_layer, 1, 1);

%% E.4.1: TFR check

im1 = q1.pgx{3};
tbaselinex = q1.tbands{1}(end-12:end-4);

for ik = 1:size(im1, 2)

    im1(:, ik, :) = im1(:, ik, :) / mean(im1(:, ik, tbaselinex), "all");

end

fmapx = q1.fmap;
locx = (linspace(1, 99, 99) - 60)*25;

figure;

tctx1 = q1.tbands{3}(end-23:end-7);
imx1 = 10*log(squeeze(mean(im1(:, :, tctx1), 3)));

subplot(2, 2, 1);
imagesc(imx1, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title("PFC (baseline before omission)");
clim([-15 15]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

tctx2 = q1.tbands{4}(1:30);
imx2 = 10*log(squeeze(mean(im1(:, :, tctx2), 3)));

subplot(2, 2, 2);
imagesc(imx2, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title("PFC (omission)");
clim([-15 15]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

subplot(2, 1, 2);
imagesc(imx2 - imx1, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title("PFC (omission - pre-omission-base)");
clim([-15 15]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

%% E4.1: Save object

temp_filename = char(q1.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q1.areainf;
q1.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q1 = load("OGLOobj\sub-C31o_ses-230816PFC.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230816PFC.mat";
tfrx = q1.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q1.jTFRplot(11, 4, q1.tbands{3}(end-10:end));
% q1.jTFRplot(3, 3, q1.tbands{2}(end-6:end), [1000 3000]);
% q1.jTFRplot(3, 4, q1.tbands{2}(end-6:end), [1000 3000]);

%% E.6: PEV calculations all trials

[expvars, layerinf] = q1.jCalcPEV(4, [1, 5]);

%% E.7: Visualize PEV

q1.jPEVplot(expvars, layerinf, [1 5]);

%% E.8: PEV calculations for all omission identities (bar plot and time plot)

% TODO: Bar plotter
% TODO: X-Y 2D TFR scatter for AX, BX, RX; A-B; A-G; B-G

[a1, a2] = q1.jTFR2dScatters(3, 7, 4, q1.tbands{2}(end-10:end), [1000 3000]);

timew = [q1.tbands{1}, q1.tbands{2}, q1.tbands{3}, q1.tbands{4}, q1.tbands{5}];
x1 = a1(2, timew);
x2 = a2(2, timew);
clx = [ones(size(q1.tbands{1})), 2*ones(size(q1.tbands{2})), ...
    3*ones(size(q1.tbands{3})), 4*ones(size(q1.tbands{4})), 5*ones(size(q1.tbands{5}))];

x1 = decimate(x1, 1);
x2 = decimate(x2, 1);

% clx = linspace(q1.tmap(1), q1.tmap(end), length(x1));

scatter(x1, x2, clx*10, clx, "filled");
colorbar;

% TODO: Save for each

% TODO: 10PFCmain.m file (and other areas)

% q1.jcalcPEVs();

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 4;
condinf = [2, 3, 4];

layerinf2 = q1.layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = q1.pgx{condinf(1), layerid}(:, q1.fbands{fband}, [q1.tbands{2}, q1.tbands{3}]); % S2d2
    x2 = q1.pgx{condinf(2), layerid}(:, q1.fbands{fband}, [q1.tbands{3}, q1.tbands{4}]); % S3d3
    x3 = q1.pgx{condinf(3), layerid}(:, q1.fbands{fband}, [q1.tbands{4}, q1.tbands{5}]); % S4d4
    
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

layerid = 4;
condinf = [6, 7, 8];

layerinf3 = q1.layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = q1.pgx{condinf(1), layerid}(:, q1.fbands{fband}, [q1.tbands{2}, q1.tbands{3}]); % S2d2
    x2 = q1.pgx{condinf(2), layerid}(:, q1.fbands{fband}, [q1.tbands{3}, q1.tbands{4}]); % S3d3
    x3 = q1.pgx{condinf(3), layerid}(:, q1.fbands{fband}, [q1.tbands{4}, q1.tbands{5}]); % S4d4
    
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

layerid = 4;
condinf = [10, 11, 12];

layerinf4 = q1.layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = q1.pgx{condinf(1), layerid}(:, q1.fbands{fband}, [q1.tbands{2}, q1.tbands{3}]); % S2d2
    x2 = q1.pgx{condinf(2), layerid}(:, q1.fbands{fband}, [q1.tbands{3}, q1.tbands{4}]); % S3d3
    x3 = q1.pgx{condinf(3), layerid}(:, q1.fbands{fband}, [q1.tbands{4}, q1.tbands{5}]); % S4d4
    
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
ntmap = q1.tmap([q1.tbands{3}, q1.tbands{4}]);

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
    xline(2000);
    xlim([min(ntmap), max(ntmap)]);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(q1.fbandlabels(fband));

end

sgtitle("Area:" + q1.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q1.freqres) + "Hz/ovlrp=." + num2str(q1.overlap) + " " + layerinf2);

%% E.13: PEV plot BX

figure;
ntmap = q1.tmap([q1.tbands{3}, q1.tbands{4}]);

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
    xline(2000);
    xlim([min(ntmap), max(ntmap)]);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(q1.fbandlabels(fband));

end

sgtitle("Area:" + q1.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q1.freqres) + "Hz/ovlrp=." + num2str(q1.overlap) + " " + layerinf3);

%% E.14: PEV plot RX

figure;
ntmap = q1.tmap([q1.tbands{3}, q1.tbands{4}]);

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
    xline(2000);
    xlim([min(ntmap), max(ntmap)]);
    xlabel("Time(ms)");ylabel("PEV(%)");
    title(q1.fbandlabels(fband));

end

sgtitle("Area:" + q1.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q1.freqres) + "Hz/ovlrp=." + num2str(q1.overlap) + " " + layerinf4);

%% Probe B (V4/MT, It seems to be majorly MT due to strong gamma response)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{2};

% q2 = load("OGLOobj\sub-C31o_ses-230816V4-MT.mat", "obj").obj;

%% E.0.1: jNWB object

q2 = jnwb(nwbFile, "V4-MT/", 500, 4250, 1, 1);

%% E.0.2: MUA plot

% q2.jMUAplot(11, [-400 4200]);
q2.jMUAplot(12, [2500 4500]);
% q2.jMUAplotAll([1000 3000]);

%% E.0.3: SUA plot

% q2.jSUAplot(3, [1000 3000], 109:110, 1:1);
q2.jSUAplot(3, [1000 1500], 109:110, 1:1);

%% E.0.4: Rastrogram

q2.jRastrogram(3, [1000 3000], 10:120, 1:2);

%% E.1: Channel and layer specs

channel_in_layer = struct();
channel_in_layer.deep = 1:27;
channel_in_layer.mid = 29:2:35;
channel_in_layer.sup = [36:44, 46:60];
channel_in_layer.goodch = [channel_in_layer.deep, channel_in_layer.mid, channel_in_layer.sup];

channel_in_layer2 = struct();
channel_in_layer2.deep = 81:128;
channel_in_layer2.mid = 76:80;
channel_in_layer2.sup = [57:64, 66:75];
channel_in_layer2.goodch = [channel_in_layer2.sup, channel_in_layer2.mid, channel_in_layer2.deep];

q2.channelinfo{1} = channel_in_layer;
q2.channelinfo{2} = channel_in_layer2;

%% E.2: LFP info plot

q2.jLFPprobeINFO(channel_in_layer.goodch);
q2.jLFPprobeINFO(channel_in_layer2.goodch);

%% E.3: Evaluate vFLIP

a1 = q2.jVFLIP(channel_in_layer.goodch, 3601:4200, 12);
a2 = q2.jVFLIP(channel_in_layer.goodch, 4101:4700, 12);

% imx1 = a1.relpow - a2.relpow;
% imx1 = smoothdata2(imx1, "movmedian", "omitmissing", "SmoothingFactor", 0.25);
% figure;
% imagesc(imx1, "XData", linspace(0, 150, size(imx1, 1)));
% clim([-.25 .25]);
% xlabel("Frequency");
% ylabel("Channel");
% title("Omission - Baseline");

%% E.3.2:

a1 = q2.jVFLIP(channel_in_layer2.goodch, 3601:4200, 12);
a2 = q2.jVFLIP(channel_in_layer2.goodch, 4101:4700, 12);

imx1 = a1.relpow - a2.relpow;
imx1 = smoothdata2(imx1, "movmedian", "omitmissing", "SmoothingFactor", 0.2);
figure;
imagesc(imx1, "XData", linspace(0, 150, size(imx1, 1)));
clim([-.25 .25]);
xlabel("Frequency");
ylabel("Channel");
title("Omission - Baseline");

%% E.4: TFR calculations all trials

% q2.jCalcTFRs(channel_in_layer, 1, 1);
q2.jCalcTFRs(channel_in_layer2, 1, 1);

%% E.4.1: TFR check

areaname = "V4";

% im1 = q2.pgx2{3} + q2.pgx2{7} + q2.pgx2{11};
im1 = q2.pgx{3} + q2.pgx{7} + q2.pgx{11};
% im1 = q2.pgx2{3};

tbaselinex = q2.tbands{1}(end-12:end-4);

for ik = 1:size(im1, 2)

    im1(:, ik, :) = im1(:, ik, :) / mean(im1(:, ik, tbaselinex), "all");

end

fmapx = q2.fmap;
locx = (linspace(1, 55, 55) - 33)*40;

figure;

tctx1 = q2.tbands{3}(end-23:end-7);
imx1 = 10*log(squeeze(mean(im1(:, :, tctx1), 3)));

subplot(2, 2, 1);
imagesc(imx1, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title(areaname + " (baseline before omission)");
clim([-15 15]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

tctx2 = q2.tbands{4}(1:30);
imx2 = 10*log(squeeze(mean(im1(:, :, tctx2), 3)));

subplot(2, 2, 2);
imagesc(imx2, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title(areaname + " (omission)");
clim([-15 15]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

subplot(2, 1, 2);
imagesc(imx2 - imx1, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title(areaname + " (omission - pre-omission-base)");
clim([-10 10]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

%% E4.1: Save object

temp_filename = char(q2.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q2.areainf;
q2.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q2 = load("OGLOobj\sub-C31o_ses-230816V4-MT.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230816V4.mat";
tfrx = q2.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

tfrpath = "tfrData\";
tfrname = "230816MT.mat";
tfrx = q2.pgx2;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q2.jTFRplot(11, 4, q2.tbands{3}(end-10:end-1));
% q2.jTFRplot(3, 1, q2.tbands{2}(end-6:end), [1000 3000]);
% q2.jTFRplot(3, 3, q2.tbands{2}(end-6:end), [1000 3000]);
% q2.jTFRplot(3, 4, q2.tbands{3}(end-10:end), [1000 3000]);
% q2.jTFRplot(3, 4, q2.tbands{2}(end-6:end), [1000 3000]);

%% E.6: PEV calculations all trials

[expvars, layerinf] = q2.jCalcPEV(1, [2, 6]);

%% E.7: Visualize PEV

q2.jPEVplot(expvars, layerinf, [2 6]);

%% E.8: PEV calculations for all omission identities (bar plot and time plot)

% q2.jcalcPEVs();

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 1;
condinf = [2, 3, 4];

layerinf2 = q2.layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = q2.pgx{condinf(1), layerid}(:, q2.fbands{fband}, q2.tbands{3}); % S2d2
    x2 = q2.pgx{condinf(2), layerid}(:, q2.fbands{fband}, q2.tbands{4}); % S3d3
    x3 = q2.pgx{condinf(3), layerid}(:, q2.fbands{fband}, q2.tbands{5}); % S4d4
    
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

layerinf3 = q2.layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = q2.pgx{condinf(1), layerid}(:, q2.fbands{fband}, q2.tbands{3}); % S2d2
    x2 = q2.pgx{condinf(2), layerid}(:, q2.fbands{fband}, q2.tbands{4}); % S3d3
    x3 = q2.pgx{condinf(3), layerid}(:, q2.fbands{fband}, q2.tbands{5}); % S4d4
    
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

layerinf4 = q2.layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = q2.pgx{condinf(1), layerid}(:, q2.fbands{fband}, q2.tbands{3}); % S2d2
    x2 = q2.pgx{condinf(2), layerid}(:, q2.fbands{fband}, q2.tbands{4}); % S3d3
    x3 = q2.pgx{condinf(3), layerid}(:, q2.fbands{fband}, q2.tbands{5}); % S4d4
    
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
ntmap = q2.tmap([q2.tbands{2}, q2.tbands{3}]);

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
    title(q2.fbandlabels(fband));

end

sgtitle("Area:" + q2.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q2.freqres) + "Hz/ovlrp=." + num2str(q2.overlap) + " " + layerinf2);

%% E.13: PEV plot BX

figure;
ntmap = q2.tmap(q2.tbands{2});

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
    title(q2.fbandlabels(fband));

end

sgtitle("Area:" + q2.areainf + " posOmission/Bx/PEV/TFR/+-2SEM/fRes=" + num2str(q2.freqres) + "Hz/ovlrp=." + num2str(q2.overlap) + " " + layerinf3);

%% E.14: PEV plot RX

figure;
ntmap = q2.tmap(q2.tbands{2});

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
    title(q2.fbandlabels(fband));

end

sgtitle("Area:" + q2.areainf + " posOmission/Rx/PEV/TFR/+-2SEM/fRes=" + num2str(q2.freqres) + "Hz/ovlrp=." + num2str(q2.overlap) + " " + layerinf4);

%% Probe C (V1/V2, It seems that it consists only of one area, to be identified, highly foveal)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{2};
% q3 = load("OGLOobj\sub-C31o_ses-230816V1-V2.mat", "obj").obj;

%% E.0.1: jNWB object

q3 = jnwb(nwbFile, "V1-V2/", 500, 4250, 2, 0);

%% E.0.2: MUA plot

q3.jMUAplot(3, [-400 4000]);

%% E.0.3: SUA plot

q3.jSUAplot(11, [100 4000], 10:20);

%% E.1: Channel and layer specs

channel_in_layer = struct(); % V1
channel_in_layer.deep = 46:107;
channel_in_layer.mid = 41:45;
channel_in_layer.sup = 11:40;
channel_in_layer.goodch = [channel_in_layer.sup, channel_in_layer.mid, channel_in_layer.deep];

q3.channelinfo{1} = channel_in_layer;

%% E.2: LFP info plot

q3.jLFPprobeINFO(channel_in_layer.goodch);

%% E.3: Evaluate vFLIP

q3.jVFLIP(channel_in_layer.goodch);

%% E.3.1: Compare vFLIP

a1 = q3.jVFLIP(channel_in_layer.goodch, 3601:4200, 12, 0.04, 0);
a2 = q3.jVFLIP(channel_in_layer.goodch, 4101:4700, 12, 0.04, 0);

imx1 = a1.relpow - a2.relpow;
imx1 = smoothdata2(imx1, "movmedian", "omitmissing", "SmoothingFactor", 0.25);
figure;
imagesc(imx1, "XData", linspace(0, 150, size(imx1, 1)));
clim([-.25 .25]);
colorbar;
xlabel("Frequency");
ylabel("Channel");
title("Omission - Baseline - V1");

%% E.4: TFR calculations all trials

% q3.jCalcTFRs(channel_in_layer);
q3.jCalcTFRs(channel_in_layer, 1, 1);

%% E.4.1: TFR check

areaname = "V1";

% im1 = q3.pgx2{3} + q3.pgx2{7} + q3.pgx2{11};
im1 = q3.pgx{3} + q3.pgx{7} + q3.pgx{11};
% im1 = q3.pgx2{3};

tbaselinex = q3.tbands{1}(end-12:end-4);

for ik = 1:size(im1, 2)

    im1(:, ik, :) = im1(:, ik, :) / mean(im1(:, ik, tbaselinex), "all");

end

fmapx = q3.fmap;
locx = (linspace(1, 55, 55) - 33)*40;

figure;

tctx1 = q3.tbands{3}(end-23:end-7);
imx1 = 10*log(squeeze(mean(im1(:, :, tctx1), 3)));

subplot(2, 2, 1);
imagesc(imx1, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title(areaname + " (baseline before omission)");
clim([-15 15]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

tctx2 = q3.tbands{4}(1:30);
imx2 = 10*log(squeeze(mean(im1(:, :, tctx2), 3)));

subplot(2, 2, 2);
imagesc(imx2, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title(areaname + " (omission)");
clim([-15 15]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

subplot(2, 1, 2);
imagesc(imx2 - imx1, "XData", fmapx, "YData", locx);
yline(0);
xlabel("Freq.");
ylabel("Dist. from L4 in um");
title(areaname + " (omission - pre-omission-base)");
clim([-10 10]);
set(gca, "YDir", "normal");
cb = colorbar();
ylabel(cb, "Power vs. baseline (dB)");

%% E4.1: Save object

temp_filename = char(q3.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q3.areainf;
q3.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q3 = load("OGLOobj\sub-C31o_ses-230816V1-V2.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230816V1.mat";
tfrx = q3.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q3.jTFRplot(3, 4, q3.tbands{3}(end-5:end), [1000 3000]);

%% E.6: PEV calculations all trials

[expvars, layerinf] = q3.jCalcPEV(1, [3, 7]);

%% E.7: Visualize PEV

q3.jPEVplot(expvars, layerinf, [3 7]);

%% E.8: PEV calculations for all omission identities (bar plot and time plot)

xs10 = q3.xs{12};
xs10a = squeeze(mean(xs10, 1));
xs10b = mean(xs10a, 1);
xs10b = smooth(xs10b, 100, "sgolay");
plot(xs10b);
% q3.jcalcPEVs();

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 1;
condinf = [2, 3, 4];

layerinf2 = q3.layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = q3.pgx{condinf(1), layerid}(:, q3.fbands{fband}, q3.tbands{3}); % S2d2
    x2 = q3.pgx{condinf(2), layerid}(:, q3.fbands{fband}, q3.tbands{4}); % S3d3
    x3 = q3.pgx{condinf(3), layerid}(:, q3.fbands{fband}, q3.tbands{5}); % S4d4
    
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

layerinf3 = q3.layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = q3.pgx{condinf(1), layerid}(:, q3.fbands{fband}, q3.tbands{3}); % S2d2
    x2 = q3.pgx{condinf(2), layerid}(:, q3.fbands{fband}, q3.tbands{4}); % S3d3
    x3 = q3.pgx{condinf(3), layerid}(:, q3.fbands{fband}, q3.tbands{5}); % S4d4
    
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

layerinf4 = q3.layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = q3.pgx{condinf(1), layerid}(:, q3.fbands{fband}, q3.tbands{3}); % S2d2
    x2 = q3.pgx{condinf(2), layerid}(:, q3.fbands{fband}, q3.tbands{4}); % S3d3
    x3 = q3.pgx{condinf(3), layerid}(:, q3.fbands{fband}, q3.tbands{5}); % S4d4
    
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
ntmap = q3.tmap(q3.tbands{2});

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
    title(q3.fbandlabels(fband));

end

sgtitle("Area:" + q3.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q3.freqres) + "Hz/ovlrp=." + num2str(q3.overlap) + " " + layerinf2);

%% E.13: PEV plot BX

figure;
ntmap = q3.tmap(q3.tbands{2});

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
    title(q3.fbandlabels(fband));

end

sgtitle("Area:" + q3.areainf + " posOmission/Bx/PEV/TFR/+-2SEM/fRes=" + num2str(q3.freqres) + "Hz/ovlrp=." + num2str(q3.overlap) + " " + layerinf3);

%% E.14: PEV plot RX

figure;
ntmap = q3.tmap(q3.tbands{2});

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
    title(q3.fbandlabels(fband));

end

sgtitle("Area:" + q3.areainf + " posOmission/Rx/PEV/TFR/+-2SEM/fRes=" + num2str(q3.freqres) + "Hz/ovlrp=." + num2str(q3.overlap) + " " + layerinf4);

%% TESTBENCH

% nwbFile = nwbPath + nwbFiles{14};
% 
% %% E.0.1: jNWB object
% 
% q2 = jnwb(nwbFile, "-/", 500, 4250, 0, 0, 0);
% 
% %% E.0.2: MUA plot
% 
% % q2.jMUAplot(11, [-400 4200]);
% q2.jMUAplotAll([-400 4000]);

%%


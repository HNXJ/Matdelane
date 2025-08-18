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

%% Ses: sub-C31o_ses-230823 (Probes A,B,C)

%% Probe A (FEF, laminar, seems like alpha/beta carry omission info positionally)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{4};

%% E.0.1: jNWB object

q1 = jnwb(nwbFile, "FEF/", 500, 4250, 0, 0);

%% E.0.2: MUA plot

q1.jMUAplot(12, [-500 4000], 100, 1);

%% E.0.3: SUA plot

q1.jSUAplot(9, [100 4000], 100:120);

%% E.1: Channel and layer specs

channel_in_layer = struct();
channel_in_layer.sup = 30:55;
channel_in_layer.mid = 56:60;
channel_in_layer.deep = 61:101;
channel_in_layer.goodch = [channel_in_layer.sup, channel_in_layer.mid, channel_in_layer.deep];

q1.channelinfo{1} = channel_in_layer;

%% E.1.1: Save single units and muas

xset = struct();
xset.xs = q1.xs;
xset.chids = q1.cs{1}.ids;
xset.peakch = q1.cs{1}.peaks;
xset.lfpch = q1.channelinfo{1};

fname = "09_FEF_convspk_1.mat";
save("spkSet\" + fname, "xset", "-v7.3");

%% E.1.2: Save LFPs

xset = struct();
xset.xs = q1.x;
xset.mdata = q1.c{1};
xset.lfpch = q1.channelinfo{1};

fname = "09_FEF_lfp_1.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

%% E.2: LFP info plot

q1.jLFPprobeINFO(channel_in_layer.goodch);
q1.jVFLIP(q1.channelinfo{1}.goodch, 3701:4200, 12);

%% E.3.0: Load Q

q1 = load("OGLOobj\sub-C31o_ses-230823FEF.mat", "obj").obj;

%% E.3: Evaluate vFLIP

a1 = q1.jVFLIP(q1.channelinfo{1}.goodch, 3601:4200, 12);
a2 = q1.jVFLIP(q1.channelinfo{1}.goodch, 4101:4700, 12);

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

q1.jCalcTFRs(channel_in_layer, 1, 1);

%%

tset = struct();
tset.pgx = q1.pgx;
tset.fmap = q1.fmap;
tset.tmap = q1.tmap;
tset.chan = q1.channelinfo;

fname = "09_FEF_tFRch_1.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

%% E4.1: Save object

temp_filename = char(q1.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q1.areainf;
q1.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q1 = load("OGLOobj\sub-C31o_ses-230823FEF.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230823FEF.mat";
tfrx = q1.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q1.jTFRplot(9, 4, q1.tbands{1}(end-10:end));
% q1.jTFRplot(10, 4, q1.tbands{1}(end-10:end));
% q1.jTFRplot(11, 4, q1.tbands{1}(end-10:end));
% q1.jTFRplot(12, 4, q1.tbands{1}(end-10:end));

%% E.6: PEV calculations all trials

[expvars, layerinf] = q1.jCalcPEV(1, [1, 5]);

%% E.7: Visualize PEV

q1.jPEVplot(expvars, layerinf, [1 5]);

%% E.8: PEV calculations for all omission identities (bar plot and time plot)

% q1.jcalcPEVs();

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 1;
condinf = [2, 3, 4];

layerinf2 = q1.layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = q1.pgx{condinf(1), layerid}(:, q1.fbands{fband}, q1.tbands{3}); % S2d2
    x2 = q1.pgx{condinf(2), layerid}(:, q1.fbands{fband}, q1.tbands{4}); % S3d3
    x3 = q1.pgx{condinf(3), layerid}(:, q1.fbands{fband}, q1.tbands{5}); % S4d4
    
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

layerinf3 = q1.layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = q1.pgx{condinf(1), layerid}(:, q1.fbands{fband}, q1.tbands{3}); % S2d2
    x2 = q1.pgx{condinf(2), layerid}(:, q1.fbands{fband}, q1.tbands{4}); % S3d3
    x3 = q1.pgx{condinf(3), layerid}(:, q1.fbands{fband}, q1.tbands{5}); % S4d4
    
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

layerinf4 = q1.layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = q1.pgx{condinf(1), layerid}(:, q1.fbands{fband}, q1.tbands{3}); % S2d2
    x2 = q1.pgx{condinf(2), layerid}(:, q1.fbands{fband}, q1.tbands{4}); % S3d3
    x3 = q1.pgx{condinf(3), layerid}(:, q1.fbands{fband}, q1.tbands{5}); % S4d4
    
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
ntmap = q1.tmap(q1.tbands{2});

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
    title(q1.fbandlabels(fband));

end

sgtitle("Area:" + q1.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q1.freqres) + "Hz/ovlrp=." + num2str(q1.overlap) + " " + layerinf2);

%% E.13: PEV plot BX

figure;
ntmap = q1.tmap(q1.tbands{2});

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
    title(q1.fbandlabels(fband));

end

sgtitle("Area:" + q1.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q1.freqres) + "Hz/ovlrp=." + num2str(q1.overlap) + " " + layerinf3);

%% E.14: PEV plot RX

figure;
ntmap = q1.tmap(q1.tbands{2});

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
    title(q1.fbandlabels(fband));

end

sgtitle("Area:" + q1.areainf + " posOmission/Ax/PEV/TFR/+-2SEM/fRes=" + num2str(q1.freqres) + "Hz/ovlrp=." + num2str(q2.overlap) + " " + layerinf4);

%% Probe B (MT/MST, It seems to be majorly MT due to strong gamma response)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{4};

%% E.0.1: jNWB object

q2 = jnwb(nwbFile, "MST/", 500, 4250, 1, 0);

%% E.0.2: MUA plot

q2.jMUAplot(2, [-400 4200]);

%% E.0.3: SUA plot

q2.jSUAplot(9, [100 4000], 100:120);

%% E.1: Channel and layer identification

channel_in_layer = struct();
channel_in_layer.deep = 73:2:110;
channel_in_layer.mid = 67:2:72;
channel_in_layer.sup = [41:2:43, 43, 47:2:64];
channel_in_layer.goodch = [channel_in_layer.sup, channel_in_layer.mid, channel_in_layer.deep];

q2.channelinfo{1} = channel_in_layer;

%% E.1.1: Save single units and muas

xset = struct();
xset.xs = q2.xs;
xset.chids = q2.cs{2}.ids;
xset.peakch = q2.cs{2}.peaks;
xset.lfpch = q2.channelinfo{1};

fname = "07_MST_convspk_2.mat";
save("spkSet\" + fname, "xset", "-v7.3");

%% E.1.2: Save LFPs

xset = struct();
xset.xs = q2.x;
xset.mdata = q2.c{1};
xset.lfpch = q2.channelinfo{1};

fname = "07_MST_lfp_2.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

%% E.2: LFP info plot

q2.jLFPprobeINFO(channel_in_layer.goodch);

%% E.3: Evaluate vFLIP

q2.jVFLIP(channel_in_layer.goodch, 1:500);

%% E.4: TFR calculations all trials

q2.jCalcTFRs(channel_in_layer, 1, 1);

%%

tset = struct();
tset.pgx = q2.pgx;
tset.fmap = q2.fmap;
tset.tmap = q2.tmap;
tset.chan = q2.channelinfo;

fname = "07_MST_tFRch_2.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

%% E4.1: Save object

temp_filename = char(q2.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q2.areainf;
q2.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q2 = load("OGLOobj\sub-C31o_ses-230823MST.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230823MST.mat";
tfrx = q2.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q2.jTFRplot(9, 4, q2.tbands{1}(end-10:end));

%% E.6: PEV calculations all trials

[expvars, layerinf] = q2.jCalcPEV(1, [2, 6]);

%% E.7: Visualize PEV

q2.jPEVplot(expvars, layerinf, [1 5]);

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
ntmap = q2.tmap(q2.tbands{2});

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

nwbFile = nwbPath + nwbFiles{4};

%% E.0.1: jNWB object

q3 = jnwb(nwbFile, "V2-V3d/", 500, 4250, 2, 0);

%% E.0.2: MUA plot

q3.jMUAplot(9, [-500 4000], 100, 1);

%% E.0.3: SUA plot

q3.jSUAplot(9, [100 4000], 100:120);

%% E.1: Channel and layer specs

channel_in_layer = struct(); % V2
channel_in_layer.deep = [45:49, 51, 53:74];
channel_in_layer.mid = 39:43;
channel_in_layer.sup = 21:38;
channel_in_layer.goodch = [channel_in_layer.sup, channel_in_layer.mid, channel_in_layer.deep];

q3.channelinfo{1} = channel_in_layer;

channel_in_layer2 = struct(); % V3d
channel_in_layer2.sup = 109:121;
channel_in_layer2.mid = 104:107;
channel_in_layer2.deep = 80:103;
channel_in_layer2.goodch = [channel_in_layer2.deep, channel_in_layer2.mid, channel_in_layer2.sup];

q3.channelinfo{2} = channel_in_layer2;

%% E.1.1: Save single units and muas

xset = struct();
xset.xs = q3.xs;
xset.chids = q3.cs{3}.ids;
xset.peakch = q3.cs{3}.peaks;
xset.lfpch = q3.channelinfo{1};

fname = "02_V2_convspk_1.mat";
save("spkSet\" + fname, "xset", "-v7.3");

xset = struct();
xset.xs = q3.xs;
xset.chids = q3.cs{3}.ids;
xset.peakch = q3.cs{3}.peaks;
xset.lfpch = q3.channelinfo{2};

fname = "03_V3d_convspk_1.mat";
save("spkSet\" + fname, "xset", "-v7.3");

%% E.1.2: Save LFPs

xset = struct();
xset.xs = q3.x;
xset.mdata = q3.c{3};
xset.lfpch = q3.channelinfo{1};

fname = "02_V2_lfp_1.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

xset = struct();
xset.xs = q3.x;
xset.mdata = q3.c{3};
xset.lfpch = q3.channelinfo{2};

fname = "03_V3d_lfp_1.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

%% E.2: LFP info plot

q3.jLFPprobeINFO(channel_in_layer.goodch);
q3.jLFPprobeINFO(channel_in_layer2.goodch);

%% E.3: Evaluate vFLIP

q3.jVFLIP(channel_in_layer.goodch);
q3.jVFLIP(channel_in_layer2.goodch);

%% E.4: TFR calculations all trials

q3.jCalcTFRs(channel_in_layer, 1, 1);
q3.jCalcTFRs(channel_in_layer2, 1, 1);

%%

tset = struct();
tset.pgx = q3.pgx;
tset.fmap = q3.fmap;
tset.tmap = q3.tmap;
tset.chan = q3.channelinfo;

fname = "02_V2_tFRch_1.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

tset = struct();
tset.pgx = q3.pgx2;
tset.fmap = q3.fmap;
tset.tmap = q3.tmap;
tset.chan = q3.channelinfo;

fname = "03_V3d_tFRch_1.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

%% E4.1: Save object

temp_filename = char(q3.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q3.areainf;
q3.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q3 = load("OGLOobj\sub-C31o_ses-230823V2-V3d.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230823V2.mat";
tfrx = q3.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

tfrpath = "tfrData\";
tfrname = "230823V3d.mat";
tfrx = q3.pgx2;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q3.jTFRplot(11, 4, q3.tbands{1}(end-10:end));

%% E.6: PEV calculations all trials

[expvars, layerinf] = q3.jCalcPEV(1, [2, 6]);

%% E.7: Visualize PEV

q3.jPEVplot(expvars, layerinf, [2 6]);

%% E.8: PEV calculations for all omission identities (bar plot and time plot)

% q3.jcalcPEVs();

%% E.9: Omission PEV (Positional, AX2/AX3/AX4)

layerid = 4;
condinf = [2, 3, 4];

layerinf2 = q3.layeridlabel(layerid) + " layer";
expvars2 = cell(1, 5);

for fband = 1:5

    x1 = q3.pgx{condinf(1), layerid}(:, q3.fbands{fband}, [q3.tbands{2}, q3.tbands{3}]); % S2d2
    x2 = q3.pgx{condinf(2), layerid}(:, q3.fbands{fband}, [q3.tbands{3}, q3.tbands{4}]); % S3d3
    x3 = q3.pgx{condinf(3), layerid}(:, q3.fbands{fband}, [q3.tbands{4}, q3.tbands{5}]); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    % expv(p > 0.05) = 1e-3;
    expv = expv .* (1 - p);
    % expv = expv .* (1 - p*2) .* (p < 0.5);
    expvars2{fband} = squeeze(expv)*100;

end

%% E.10: Omission PEV (Positional, BX2/BX3/BX4)

layerid = 4;
condinf = [6, 7, 8];

layerinf3 = q3.layeridlabel(layerid) + " layer";
expvars3 = cell(1, 5);

for fband = 1:5

    x1 = q3.pgx{condinf(1), layerid}(:, q3.fbands{fband}, [q3.tbands{2}, q3.tbands{3}]); % S2d2
    x2 = q3.pgx{condinf(2), layerid}(:, q3.fbands{fband}, [q3.tbands{3}, q3.tbands{4}]); % S3d3
    x3 = q3.pgx{condinf(3), layerid}(:, q3.fbands{fband}, [q3.tbands{4}, q3.tbands{5}]); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    % expv(p > 0.05) = 1e-3;
    expv = expv .* (1 - p);
    % expv = expv .* (1 - p*2) .* (p < 0.5);
    expvars3{fband} = squeeze(expv)*100;

end

%% E.11: Omission PEV (Positional, RX2/RX3/RX4)

layerid = 4;
condinf = [10, 11, 12];

layerinf4 = q3.layeridlabel(layerid) + " layer";
expvars4 = cell(1, 5);

for fband = 1:5

    x1 = q3.pgx{condinf(1), layerid}(:, q3.fbands{fband}, [q3.tbands{2}, q3.tbands{3}]); % S2d2
    x2 = q3.pgx{condinf(2), layerid}(:, q3.fbands{fband}, [q3.tbands{3}, q3.tbands{4}]); % S3d3
    x3 = q3.pgx{condinf(3), layerid}(:, q3.fbands{fband}, [q3.tbands{4}, q3.tbands{5}]); % S4d4
    
    [N1, nF, nT] = size(x1);
    N2 = size(x2, 1);
    N3 = size(x3, 1);
    
    data = zeros(N1+N2+N3, nF, nT);
    data(1:N1, :, :) = x1;
    data(N1+1:N1+N2, :, :) = x2;
    data(N1+N2+1:N1+N2+N3, :, :) = x3;
    
    groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    % expv(p > 0.05) = 1e-3;
    expv = expv .* (1 - p);
    % expv = expv .* (1 - p*2) .* (p < 0.5);
    expvars4{fband} = squeeze(expv)*100;

end

%% E.12: PEV plot AX

figure;
ntmap = q3.tmap([q3.tbands{3}, q3.tbands{4}]);

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
ntmap = q3.tmap([q3.tbands{3}, q3.tbands{4}]);

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
ntmap = q3.tmap([q3.tbands{3}, q3.tbands{4}]);

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
    title(q3.fbandlabels(fband));

end

sgtitle("Area:" + q3.areainf + " posOmission/Rx/PEV/TFR/+-2SEM/fRes=" + num2str(q3.freqres) + "Hz/ovlrp=." + num2str(q3.overlap) + " " + layerinf4);

%%
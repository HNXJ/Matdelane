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

%% Ses: sub-V198o_ses-230720 (Probes A,B)

%% Probe A (V1-V2, laminar)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{15};

%% E.0.1: jNWB object

q1 = jnwb(nwbFile, "V1-V2/", 500, 4250, 0, 0);

%% E.0.2: MUA plot

q1.jMUAplot(7, [-500 4000]);

%% E.0.3: SUA plot

q1.jSUAplot(9, [100 4000], 100:120);

%% E.1: Channel and layer specs

channel_in_layer = struct(); % V1
channel_in_layer.sup = 1:11;
channel_in_layer.mid = 12:15;
channel_in_layer.deep = 16:34;
channel_in_layer.goodch = [channel_in_layer.sup, channel_in_layer.mid, channel_in_layer.deep];

channel_in_layer2 = struct(); % V2
channel_in_layer2.deep = [58:65, 65, 67:74];
channel_in_layer2.mid = 75:78;
channel_in_layer2.sup = 79:95;
channel_in_layer2.goodch = [channel_in_layer2.deep, channel_in_layer2.mid, channel_in_layer2.sup];

q1.channelinfo{1} = channel_in_layer;
q1.channelinfo{2} = channel_in_layer2;

%% E.1.1: Save single units and muas

xset = struct();
xset.xs = q1.xs;
xset.chids = q1.cs{1}.ids;
xset.peakch = q1.cs{1}.peaks;
xset.lfpch = q1.channelinfo{1};

fname = "01_V1_convspk_5.mat";
save("spkSet\" + fname, "xset", "-v7.3");

xset = struct();
xset.xs = q1.xs;
xset.chids = q1.cs{1}.ids;
xset.peakch = q1.cs{1}.peaks;
xset.lfpch = q1.channelinfo{2};

fname = "02_V2_convspk_4.mat";
save("spkSet\" + fname, "xset", "-v7.3");

% E.1.2: Save LFPs

xset = struct();
xset.xs = q1.x;
xset.mdata = q1.c{1};
xset.lfpch = q1.channelinfo{1};

fname = "01_V1_lfp_5.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

xset = struct();
xset.xs = q1.x;
xset.mdata = q1.c{2};
xset.lfpch = q1.channelinfo{2};

fname = "02_V2_lfp_4.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

%% E.2: LFP info plot

q1.jLFPprobeINFO(channel_in_layer.goodch);
q1.jLFPprobeINFO(channel_in_layer2.goodch);

%% E.3: Evaluate vFLIP

q1.jVFLIP(channel_in_layer.goodch);
q1.jVFLIP(channel_in_layer2.goodch);

%% E.4: TFR calculations all trials

q1.jCalcTFRs(channel_in_layer, 1, 1);
q1.jCalcTFRs(channel_in_layer2, 1, 1);

%%

tset = struct();
tset.pgx = q1.pgx;
tset.fmap = q1.fmap;
tset.tmap = q1.tmap;
tset.chan = q1.channelinfo;

fname = "01_V1_tFRch_5.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

tset = struct();
tset.pgx = q1.pgx2;
tset.fmap = q1.fmap;
tset.tmap = q1.tmap;
tset.chan = q1.channelinfo;

fname = "02_V2_tFRch_4.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

%% E4.1: Save object

temp_filename = char(q1.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q1.areainf;
q1.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q1 = load("OGLOobj\sub-V198o_ses-230720V1-V2.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230720V1.mat";
tfrx = q1.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

tfrpath = "tfrData\";
tfrname = "230720V2.mat";
tfrx = q1.pgx2;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q1.jTFRplot(9, 4, xinfo.tbands{1}(end-5:end));

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

%% Probe B (V3d/V3a)

%% E.0: Load NWB

nwbFile = nwbPath + nwbFiles{15};

%% E.0.1: jNWB object

q2 = jnwb(nwbFile, "V3d-V3a/", 500, 4250, 1, 0);

%% E.0.2: MUA plot

q2.jMUAplot(9, [1000 3000]);

%% E.0.3: SUA plot

q2.jSUAplot(9, [100 4000], 100:120);

%% E.1: Channel and layer identification

channel_in_layer = struct(); % V3d
channel_in_layer.deep = [20:27, 27, 29, 29, 31:40];
channel_in_layer.mid = 16:19;
channel_in_layer.sup = [8:9, 9, 11, 11, 13, 13, 15, 15];
channel_in_layer.goodch = [channel_in_layer.sup, channel_in_layer.mid, channel_in_layer.deep];

channel_in_layer2 = struct(); % V3a
channel_in_layer2.deep = 1 + (41:74) - mod(41:74, 2);
channel_in_layer2.mid = 1 + (75:79) - mod(75:79, 2);
channel_in_layer2.sup = 1 + (81:95) - mod(81:95, 2);
channel_in_layer2.goodch = [channel_in_layer2.deep, channel_in_layer2.mid, channel_in_layer2.sup];

q2.channelinfo{1} = channel_in_layer;
q2.channelinfo{2} = channel_in_layer2;

%% E.1.1: Save single units and muas

xset = struct();
xset.xs = q2.xs;
xset.chids = q2.cs{2}.ids;
xset.peakch = q2.cs{2}.peaks;
xset.lfpch = q2.channelinfo{1};

fname = "03_V3d_convspk_5.mat";
save("spkSet\" + fname, "xset", "-v7.3");

xset = struct();
xset.xs = q2.xs;
xset.chids = q2.cs{2}.ids;
xset.peakch = q2.cs{2}.peaks;
xset.lfpch = q2.channelinfo{2};

fname = "04_V3a_convspk_2.mat";
save("spkSet\" + fname, "xset", "-v7.3");

% E.1.2: Save LFPs

xset = struct();
xset.xs = q2.x;
xset.mdata = q2.c{1};
xset.lfpch = q2.channelinfo{1};

fname = "03_V3d_lfp_5.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

xset = struct();
xset.xs = q2.x;
xset.mdata = q2.c{2};
xset.lfpch = q2.channelinfo{2};

fname = "04_V3a_lfp_2.mat";
save("lfpSet\" + fname, "xset", "-v7.3");

%% E.2: LFP info plot

q2.jLFPprobeINFO(channel_in_layer.goodch);
q2.jLFPprobeINFO(channel_in_layer2.goodch);

%% E.3: Evaluate vFLIP

q2.jVFLIP(channel_in_layer.goodch, 1:1000);
q2.jVFLIP(channel_in_layer2.goodch, 1:1000);

%% E.4: TFR calculations all trials

q2.jCalcTFRs(channel_in_layer, 1, 1);
q2.jCalcTFRs(channel_in_layer2, 1, 1);

%%

tset = struct();
tset.pgx = q2.pgx;
tset.fmap = q2.fmap;
tset.tmap = q2.tmap;
tset.chan = q2.channelinfo;

fname = "03_V3d_tFRch_5.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

tset = struct();
tset.pgx = q2.pgx2;
tset.fmap = q2.fmap;
tset.tmap = q2.tmap;
tset.chan = q2.channelinfo;

fname = "04_V3a_tFRch_2.mat";
save("tfrSet\" + fname, "tset", "-v7.3");

%% E4.1: Save object

temp_filename = char(q2.nwbFile);
temp_filename = temp_filename(6:end-4);
temp_filename = temp_filename + q2.areainf;
q2.jSave("OGLOobj", temp_filename);

%% E4.2: Load if object exists

q2 = load("OGLOobj\sub-V198o_ses-230720V3d-V3a.mat", "obj").obj;

%% E4.3: Save TFR separately

tfrpath = "tfrData\";
tfrname = "230720V3d.mat";
tfrx = q2.pgx;
save(tfrpath + tfrname, "tfrx", "-v7.3");

tfrpath = "tfrData\";
tfrname = "230720V3a.mat";
tfrx = q2.pgx2;
save(tfrpath + tfrname, "tfrx", "-v7.3");

%% E.5: Visualize TFR

q2.jTFRplot(9, 4, q2.tbands{1}(end-10:end));

%% E.6: PEV calculations all trials

[expvars, layerinf] = q2.jCalcPEV(1, [1, 5]);

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

%%
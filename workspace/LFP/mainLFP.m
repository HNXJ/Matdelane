%% Setup

clear;clc;close all;

cd('D:\Electrophysiology\Matdelane\');
addpath(genpath('matnwb'));
addpath(genpath('matdelane'));
addpath(genpath('flipv2'));
generateCore();

addpath('fieldtrip');ft_defaults;

disp("Toolbox setup done.");

%% E.1: Load LFP

nPath = "lfpData\";
nFiles = {dir(nPath).name};
nFiles = nFiles(endsWith(nFiles, ".mat"));
xL = cell(1, 10);

for ik = 1:10

    nFile = nFiles{ik};
    tic;
    xL{ik} = load(nPath + nFile).z;
    toc;

end

%% E.2: TFR from LFP; Theta / Alpha / Beta / GammaL / GammaH

frqlims = [0 250];
timeres = 0.2;
overlap = 95;

[p1, ~, ~] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", frqlims, "TimeResolution", timeres, "OverlapPercent", overlap);

%%

g1 = xL{9}{1, 1}{1};
g2 = xL{10}{1, 1}{6};

imx1 = squeeze(mean(g1, 2));
imx2 = squeeze(mean(g2, 2));

figure;
subplot(2, 1, 1);
imagesc(imx1);

subplot(2, 1, 2);
imagesc(imx2);

%%

[N1, N2, N3, N4] = size(xL{1}{1, 1}{1});

for ik = 1 % N1

    for jk = 1:N2

        for kk = 1:N3



        end

    end

end

%% G.1: PEV

cond1 = 1; % AAAB
cond2 = 5; % BBBA

expvars2 = cell(1, 10);

for k = 1:10

    ntotal = 0;

    for ses = 1:size(xL{k}{1, cond1}, 2)

        nchan = size(xL{k}{1, cond1}{1}, 2);
        expvars2{k} = zeros(nchan, 5000);
        disp(k);
    
        x1 = xL{k}{1, cond1}{ses}(:, :, 1:5000);
        x2 = xL{k}{1, cond2}{ses}(:, :, 1:5000);
        
        N = min(size(x1, 1), size(x2, 1));
        n1 = min(size(x1, 2), size(x2, 2));
    
        idx1 = randperm(size(x1, 1), N);
        idx2 = randperm(size(x2, 1), N);
        
        % x1 = squeeze(mean(x1(idx1, :, :), 1));
        % x2 = squeeze(mean(x2(idx2, :, :), 1));
        idxch = mod(randperm(nchan, nchan), n1) + 1;
        x1 = x1(idx1, idxch, :);
        x2 = x2(idx2, idxch, :);
        
        if n1 > 0 && N > 0 
    
            data = zeros(N*2, nchan, 5000);
            data(1:N, :, :) = x1;
            data(N+1:end, :, :) = x2;
            data = jSmooth(data, 100);
    
            groupIDs = [ones(1, N), ones(1, N)*2];
            [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
            expvars2{k}(ntotal+1:ntotal+nchan, :) = squeeze(expv);
    
        end
    
        ntotal = n1 + ntotal;

    end

end


%% G.2: Plot PEV 1-5

areas = ["V1", "V2", "V3d", "V3a", "V4", "TEO", "MT", "MST", "FEF", "PFC"];
figure("Position", [0 0 2000 1000]);
cl = [0.9 0.7 0.7];

for k = 1:5

    subplot(5, 1, mod(k-1, 5)+1);
    ths = median(mean(expvars{k}, 2));
    idxp = mean(expvars{k}, 2) > ths;
    % figure("Position", [0 0 1500 500]);
    t = linspace(-970, 4030, 5000);
    yt = detrend(squeeze(mean(expvars{k}(idxp, :), 1)))'*100;
    yt(yt<0) = 0;
    st = smooth(squeeze(std(expvars{k}(idxp, :), 1)), 50)*100/sqrt(64);
    stx = smooth(yt + st, 20);stx(stx<0) = 0;
    sty = smooth(yt - st, 20);sty(sty<0) = 0;
    plot(t, stx, "Color", cl);
    hold("on");
    plot(t, sty, "Color", cl);
    patch([t, t(end:-1:1)], [stx;sty(end:-1:1)], cl);
    plot(t, yt, "Color", "k");

    xline(-500, "HandleVisibility", "off", "Color", "g");
    xline(0, "HandleVisibility", "off", "Color", "r");
    xline(1030, "HandleVisibility", "off", "Color", "r");
    xline(2060, "HandleVisibility", "off", "Color", "r");
    xline(3090, "HandleVisibility", "off", "Color", "r");

    title(areas(k));
    xlabel("Time (ms)");
    ylabel("PEV (%)");
    xlim([-1000 4000]);

end

sgtitle("PEV; LFP; average by area; AAAB vs BBBA");
fname = num2str(k) + "_abpev";
print(gcf,'-vector','-dsvg',fname +".svg");

%% G.2: Plot PEV 6-10

areas = ["V1", "V2", "V3d", "V3a", "V4", "TEO", "MT", "MST", "FEF", "PFC"];
figure("Position", [0 0 2000 1000]);
cl = [0.9 0.7 0.7];

for k = 6:10

    subplot(5, 1, mod(k-1, 5)+1);
    ths = median(mean(expvars{k}, 2));
    idxp = mean(expvars{k}, 2) > ths;
    % figure("Position", [0 0 1500 500]);
    t = linspace(-970, 4030, 5000);
    yt = detrend(squeeze(mean(expvars{k}(idxp, :), 1)))'*100;
    yt(yt<0) = 0;
    st = smooth(squeeze(std(expvars{k}(idxp, :), 1)), 50)*100/sqrt(64);
    stx = smooth(yt + st, 20);stx(stx<0) = 0;
    sty = smooth(yt - st, 20);sty(sty<0) = 0;
    plot(t, stx, "Color", cl);
    hold("on");
    plot(t, sty, "Color", cl);
    patch([t, t(end:-1:1)], [stx;sty(end:-1:1)], cl);
    plot(t, yt, "Color", "k");

    xline(-500, "HandleVisibility", "off", "Color", "g");
    xline(0, "HandleVisibility", "off", "Color", "r");
    xline(1030, "HandleVisibility", "off", "Color", "r");
    xline(2060, "HandleVisibility", "off", "Color", "r");
    xline(3090, "HandleVisibility", "off", "Color", "r");

    title(areas(k));
    xlabel("Time (ms)");
    ylabel("PEV (%)");
    xlim([-1000 4000]);

end

sgtitle("PEV; LFP; average by area; AAAB vs BBBA");
fname = num2str(k) + "_abpev";
print(gcf,'-vector','-dsvg',fname +".svg");

%%

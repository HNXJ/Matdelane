%% Setup

clear;clc;close all;

cd('D:\Electrophysiology\Matdelane\');
addpath(genpath('matnwb'));
addpath(genpath('matdelane'));
addpath(genpath('flipv2'));
generateCore();

addpath('fieldtrip');ft_defaults;

disp("Toolbox setup done.");

%% E.1: Load TFR

nPath = "lfpData\";
nFiles = {dir(nPath).name};
nFiles = nFiles(endsWith(nFiles, ".mat"));
xG = cell(1, 3);

tic;

for ik = 16:18

    nFile = nFiles{ik};
    xG{ik-15} = load(nPath + nFile);

end

toc;

%% E.2: TODO** Compress TFR across sessions for PEV

pgxx = xG{1}.pgx;

%%



%% G.1: PEV

cond1 = 1; % AAAB
cond2 = 5; % AAAX

% expvars2 = cell(1, 10);
fcnt = size(pgxx, 4);
tcnt = size(pgxx, 5);
fmap = linspace(0, 250, fcnt);
tmap = linspace(-900, 5100, tcnt);

thetamap = find(fmap > 2, 1):find(fmap > 7, 1);
alphamap = find(fmap > 8, 1):find(fmap > 12, 1);
betamap = find(fmap > 14, 1):find(fmap > 30, 1);
gammamap1 = find(fmap > 32, 1):find(fmap > 64, 1);
gammamap2 = find(fmap > 64, 1):find(fmap >= 249, 1);
    
for k = 1 % area

    ntotal = 0;
    pgxx = xG{k}.pgx;
    bands = gammamap2;

    for ses = 1:size(pgxx, 1) % ses

        nchan = size(pgxx, 3);
        nF = size(bands, 2);
        nT = 550;
        expvars2{k} = zeros(nchan, nF, nT);
        disp(k);
    
        x1 = squeeze(pgxx(ses, cond1, :, bands, 1:nT));
        x2 = squeeze(pgxx(ses, cond2, :, bands, 1:nT));
        % x2 = squeeze(pgxx(ses, cond1, :, bands, randperm(nT, nT)));
        % imgx = squeeze(mean(mean(x1, 2), 3));
        % aimgx = imgx > max(imgx)*0.0;
        % 
        % x1(aimgx, :, :) = 0;
        % x2(aimgx, :, :) = 0;
        
        N = min(size(x1, 1), size(x2, 1));
        n1 = min(size(x1, 2), size(x2, 2));
    
        idx1 = randperm(size(x1, 1), N);
        idx2 = randperm(size(x2, 1), N);
        
        % x1 = squeeze(mean(x1(idx1, :, :), 1));
        % x2 = squeeze(mean(x2(idx2, :, :), 1));
        % idxch = mod(randperm(nchan, nchan), n1) + 1;
        % x1 = x1(idxch, :);
        % x2 = x2(idxch, :);
        
        if n1 > 0 && N > 0 
    
            data = zeros(N*2, n1, nT);
            data(1:N, :, :) = x1;
            data(N+1:end, :, :) = x2;
            data = jSmooth(data, 10);
    
            groupIDs = [ones(1, N), ones(1, N)*2];
            [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
            expvars2{k}(ntotal+1:ntotal+1, :, :) = squeeze(expv);
    
        end
    
        ntotal = n1 + ntotal;

    end

end


%% G.2: Plot PEV 1-5

areas = ["V1", "V2", "V3d", "V3a", "V4", "TEO", "MT", "MST", "FEF", "PFC"];
figure("Position", [0 0 2000 1000]);
cl = [0.9 0.7 0.7];
expvars = expvars2;

for k = 1

    subplot(3, 1, mod(k-1, 5)+1);
    expvartemp = squeeze(mean(expvars{k}, 2));
    ths = median(mean(expvartemp, 2));
    idxp = mean(expvartemp, 2) > 0;
    t = tmap(1:550);
    yt = detrend(squeeze(mean(expvartemp(idxp, :), 1)))'*100;
    yt(yt<0) = 0;
    st = smooth(squeeze(std(expvartemp(idxp, :), 1)), 10)*100/sqrt(64);
    stx = smooth(yt + st, 2);stx(stx<0) = 0;
    sty = smooth(yt - st, 2);sty(sty<0) = 0;
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

    title(areas(k+5));
    xlabel("Time (ms)");
    ylabel("PEV (%)");
    % xlim([-1000 5000]);

end

sgtitle("PEV; TFR; average by band; AAAB vs null");
fname = num2str(k) + "_abpevB";
print(gcf,'-vector','-dsvg',fname +".svg");


%% PEVs on TFRs
% Calculate PEV on Bands 
% Calculate PEV on VHG (very high gamma >150)
% Should be similar to ~sMUA or ~eMUA

%%

img = squeeze(mean(xG{5}.pgx(3, 9, 1:2:128, :, 1:500), 3, "omitnan"));
img2 = squeeze(mean(xG{5}.pgx(3, 12, 1:2:128, :, 1:500), 3, "omitnan"));
fcnt = size(img, 1);
tcnt = size(img, 2);

for ij = 1:fcnt

    img(ij, :) = smooth(img(ij, :), 40, "lowess");
    img2(ij, :) = smooth(img2(ij, :), 40, "lowess");
    img(ij, :) = img(ij, :) / mean(img(ij, 1:40));
    img2(ij, :) = img2(ij, :) / mean(img2(ij, 1:40));

    % img(ij, :) = (img(ij, :) - mean(img(ij, :))) / (std(img(ij, :))/sqrt(fcnt));
    % img2(ij, :) = (img2(ij, :) - mean(img2(ij, :))) / (std(img2(ij, :))/sqrt(fcnt));

end

fmap = linspace(0, 250, fcnt);
tmap = linspace(-900, 5100, tcnt);

subplot(2, 1, 1);
imagesc(img, "XData", tmap, "YData", fmap);
ylim([0 100]);
clim([0 5]);
set(gca, "YDir", "normal");
colorbar;

subplot(2, 1, 2);
imagesc(img2, "XData", tmap, "YData", fmap);
ylim([0 100]);
clim([0 5]);
set(gca, "YDir", "normal");
colorbar;

%%

areas = ["V1", "V2", "V3d", "V3a", "V4", "TEO", "MT", "MST", "FEF", "PFC"];

ik = 5;
jZscoreTFR(xG{ik}.pgx, areas(ik), 10, 0.1);

%%

for ik = 5

    jPlotTFR(xG{ik}.pgx, areas(ik), 10, 0.1);
    jPlotTFRhq(xG{ik}.pgx, areas(ik), 10, 0.1, 1:3);

end

%%

jPlotTFRhq(xG{8}.pgx, xG{8}.pgxcnt, areas(8), 10, 0.07, 1);
% jPlotTFRhq(xG{8}.pgx, xG{8}.pgxcnt, areas(8), 10, 0.1, 2);
% jPlotTFRhq(xG{8}.pgx, xG{8}.pgxcnt, areas(8), 50, 0.07, 3);
% jPlotTFRhq(xG{8}.pgx, xG{8}.pgxcnt, areas(8), 10, 0.1, 0);

% jPlotTFR(xG{1}.pgx, "V1", 10, 0.15);
% jPlotTFR(xG{10}.pgx, "PFC", 10, 0.15);

%% F.1: TFR plotter by condition

% TODO : 1. Select only good channels
% TODO : 2. Foveal/pFoveal/nFoveal
% TODO : 3. vFLIP pipeline

function jPlotTFR(pgx, areaName, smoothW, baselineW)
    
    fcnt = size(pgx, 4);
    tcnt = size(pgx, 5);
    fmap = linspace(0, 250, fcnt);
    tmap = linspace(-900, 5100, tcnt);

    thetamap = find(fmap > 2, 1):find(fmap > 7, 1);
    alphamap = find(fmap > 8, 1):find(fmap > 12, 1);
    betamap = find(fmap > 14, 1):find(fmap > 30, 1);
    gammamap1 = find(fmap > 32, 1):find(fmap > 64, 1);
    gammamap2 = find(fmap > 64, 1):find(fmap >= 249, 1);
    
    pgx2 = squeeze(mean(pgx, 1));
    conditionNames = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA", "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

    for ib = [0, 4, 8]

        figure;

        for ic = 1:4
    
            conditionID = ic + ib;
            conditionName = conditionNames(conditionID);
    
            subplot(4, 1, ic);
            y1t = smooth(squeeze(mean(mean(pgx2(conditionID, :, thetamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1a = smooth(squeeze(mean(mean(pgx2(conditionID, :, alphamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1b = smooth(squeeze(mean(mean(pgx2(conditionID, :, betamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1g1 = smooth(squeeze(mean(mean(pgx2(conditionID, :, gammamap1, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1g2 = smooth(squeeze(mean(mean(pgx2(conditionID, :, gammamap2, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            
            y1t = baselineCorrect(y1t, baselineW);
            y1a = baselineCorrect(y1a, baselineW);
            y1b = baselineCorrect(y1b, baselineW);
            y1g1 = baselineCorrect(y1g1, baselineW);
            y1g2 = baselineCorrect(y1g2, baselineW);
    
            plot(tmap, 10*log10((y1t)), "DisplayName", "RX-Theta");hold("on");
            plot(tmap, 10*log10((y1a)), "DisplayName", "RX-Alpha");
            plot(tmap, 10*log10((y1b)), "DisplayName", "RX-Beta");
            plot(tmap, 10*log10((y1g1)), "DisplayName", "RX-Gamma1");
            plot(tmap, 10*log10((y1g2)), "DisplayName", "RX-Gamma2");
            xline(-471, "LineStyle", "-", "HandleVisibility", "off");
            xline(-31, "LineStyle", "--", "HandleVisibility", "off");
            % xline(4.500, "LineStyle", "--", "HandleVisibility", "off");
            legend();title(areaName + "-" + conditionName);
            ylim([-15 25]);
            xlim([-900 5100]);
        end

        fname = areaName + "_" + conditionName(1) + "_TFRBandplot";
        print(gcf,'-vector','-dsvg',fname +".svg");

    end
    
end


function jPlotTFRhq(pgx, pgxcnt, areaName, smoothW, baselineW, chanx)
    
    fcnt = size(pgx, 4);
    tcnt = size(pgx, 5);
    fmap = linspace(0, 250, fcnt);
    tmap = linspace(-900, 5100, tcnt);

    thetamap = find(fmap > 2, 1):find(fmap > 7, 1);
    alphamap = find(fmap > 8, 1):find(fmap > 12, 1);
    betamap = find(fmap > 14, 1):find(fmap > 30, 1);
    gammamap1 = find(fmap > 32, 1):find(fmap > 64, 1);
    gammamap2 = find(fmap > 64, 1):find(fmap >= 249, 1);
    
    pgxc2 = zeros([size(pgx, 1), size(pgx, 2), 2, size(pgx, 4), size(pgx, 5)]);

    for ik = 1:size(pgx, 1)

        for jk = 1:12

            chanset = squeeze(pgxcnt(ik, 1, :, 2) == chanx);
            pgxc2(ik, jk, 1, :, :) = mean(pgx(ik, jk, chanset, :, :), 3, "omitnan");
            pgxc2(ik, jk, 2, :, :) = std(mean(mean(pgx(ik, jk, chanset, :, :), 4, "omitnan"), 3, "omitnan"));

        end

    end

    pgx2 = squeeze(mean(pgxc2, 1, "omitnan"));
    conditionNames = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA", "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

    for ib = [0, 4, 8]

        figure;

        for ic = 1:4
    
            conditionID = ic + ib;
            conditionName = conditionNames(conditionID);
    
            subplot(4, 1, ic);
            y1t = smooth(squeeze(mean(mean(pgx2(conditionID, 1, thetamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'lowess');
            y1a = smooth(squeeze(mean(mean(pgx2(conditionID, 1, alphamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'lowess');
            y1b = smooth(squeeze(mean(mean(pgx2(conditionID, 1, betamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'lowess');
            y1g1 = smooth(squeeze(mean(mean(pgx2(conditionID, 1, gammamap1, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'lowess');
            y1g2 = smooth(squeeze(mean(mean(pgx2(conditionID, 1, gammamap2, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'lowess');
            
            y1t = baselineCorrect(y1t, baselineW);
            y1a = baselineCorrect(y1a, baselineW);
            y1b = baselineCorrect(y1b, baselineW);
            y1g1 = baselineCorrect(y1g1, baselineW);
            y1g2 = baselineCorrect(y1g2, baselineW);
    
            plot(tmap, 10*log10((y1t)), "DisplayName", "RX-Theta");hold("on");
            plot(tmap, 10*log10((y1a)), "DisplayName", "RX-Alpha");
            plot(tmap, 10*log10((y1b)), "DisplayName", "RX-Beta");
            plot(tmap, 10*log10((y1g1)), "DisplayName", "RX-Gamma1");
            plot(tmap, 10*log10((y1g2)), "DisplayName", "RX-Gamma2");
            xline(-471, "LineStyle", "-", "HandleVisibility", "off");
            xline(-31, "LineStyle", "--", "HandleVisibility", "off");
            % xline(4.500, "LineStyle", "--", "HandleVisibility", "off");
            legend();title(areaName + "-" + conditionName);
            ylim([-10 25]);
            xlim([-500 5100]);
        end

        % fname = areaName + "_" + conditionName(1) + "_TFRBandplot";
        % print(gcf,'-vector','-dsvg',fname +".svg");

    end
    
end


function jZscoreTFR(pgx, areaName, smoothW, baselineW)
    
    fcnt = size(pgx, 4);
    tcnt = size(pgx, 5);
    fmap = linspace(0, 250, fcnt);
    tmap = linspace(-900, 5100, tcnt);

    thetamap = find(fmap > 2, 1):find(fmap > 7, 1);
    alphamap = find(fmap > 8, 1):find(fmap > 12, 1);
    betamap = find(fmap > 14, 1):find(fmap > 30, 1);
    gammamap1 = find(fmap > 32, 1):find(fmap > 64, 1);
    gammamap2 = find(fmap > 64, 1):find(fmap >= 249, 1);
    
    pgx2 = squeeze(mean(pgx, 1));
    conditionNames = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA", "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

    for ib = [0, 4, 8]

        figure;

        for ic = 1:4
    
            conditionID = ic + ib;
            conditionName = conditionNames(conditionID);
    
            subplot(4, 1, ic);
            y1t = smooth(squeeze(mean(mean(pgx2(conditionID, :, thetamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1a = smooth(squeeze(mean(mean(pgx2(conditionID, :, alphamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1b = smooth(squeeze(mean(mean(pgx2(conditionID, :, betamap, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1g1 = smooth(squeeze(mean(mean(pgx2(conditionID, :, gammamap1, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            y1g2 = smooth(squeeze(mean(mean(pgx2(conditionID, :, gammamap2, :), 2, "omitnan"), 3, "omitnan")), smoothW, 'loess');
            
            y1t = baselineCorrect(y1t, baselineW);
            y1a = baselineCorrect(y1a, baselineW);
            y1b = baselineCorrect(y1b, baselineW);
            y1g1 = baselineCorrect(y1g1, baselineW);
            y1g2 = baselineCorrect(y1g2, baselineW);
    
            plot(tmap, 10*log10((y1t)), "DisplayName", "RX-Theta");hold("on");
            plot(tmap, 10*log10((y1a)), "DisplayName", "RX-Alpha");
            plot(tmap, 10*log10((y1b)), "DisplayName", "RX-Beta");
            plot(tmap, 10*log10((y1g1)), "DisplayName", "RX-Gamma1");
            plot(tmap, 10*log10((y1g2)), "DisplayName", "RX-Gamma2");
            xline(-471, "LineStyle", "-", "HandleVisibility", "off");
            xline(-31, "LineStyle", "--", "HandleVisibility", "off");
            % xline(4.500, "LineStyle", "--", "HandleVisibility", "off");
            legend();title(areaName + "-" + conditionName);
            ylim([-15 25]);
            xlim([-900 5100]);
        end

        fname = areaName + "_" + conditionName(1) + "_TFRBandplot";
        print(gcf,'-vector','-dsvg',fname +".svg");

    end
    
end


function y = baselineCorrect(x, ratioW)

    k = ceil(ratioW * length(x));
    baselineW = mean(x(1:k));
    y = x / baselineW;

end

%% Z.1: Load LFP's and save TFRs by session and RF stats

nwbPath = "data\";
nwbFiles = {dir(nwbPath).name};
nwbFiles = nwbFiles(endsWith(nwbFiles, ".nwb"));

nPath = "lfpData\";
nFiles = {dir(nPath).name};
nFiles = nFiles(endsWith(nFiles, ".mat"));
areas = ["V1", "V2", "V3d", "V3a", "V4", "TEO", "MT", "MST", "FEF", "PFC"];

for i = 9:10

    nFile = nFiles{i};
    xG = load(nPath + nFile).z;
    disp(num2str(i) + ". loaded : " + nFile);
    [pgx, pgxcnt] = jOGLOTFRs(xG);

    fname = jFileNameTFR(i, areas(i));
    save(fname, "pgx", "pgxcnt", "-v7.3");
    clear pgx pgxcnt xG

end

%% XF.1: TFR filename maker

function fname = jFileNameTFR(acnt, area)

   if acnt < 10
    
        fname = "lfpData\TFRs_0" + num2str(acnt) + "_" + area;

    else

        fname = "lfpData\TFRs_" + num2str(acnt) + "_" + area;

    end 

end

%% EoF
%% Setup

clear;clc;close all;

cd('D:\Electrophysiology\Matdelane\');
addpath(genpath('matnwb'));
addpath(genpath('matdelane'));
addpath(genpath('flipv2'));

generateCore();
addpath('fieldtrip');
ft_defaults;
disp("Setup done.");

load("tfrSet\info.mat");

%% SPK unifier

areax = "conv";
spkpath = "spkSet\";
spkfiles = {dir(spkpath).name};
spkfiles = spkfiles(contains(spkfiles, areax));

Nfiles = length(spkfiles);
spkArea = cell(Nfiles, 1);
sspkData = cell(Nfiles, 12);

%% LFP unifier

areax = "FEF";
lfppath = "lfpSet\";
lfpfiles = {dir(lfppath).name};
lfpfiles = lfpfiles(contains(lfpfiles, areax));
Nfiles = length(lfpfiles);
lfpData = cell(Nfiles, 1);

%% Load spk

spkData = cell(Nfiles, 1);

for ik = 1:Nfiles

    tsetx = load(spkpath + spkfiles{ik});
    tsetf = strsplit(spkfiles{ik}, "_");
    spkArea{ik} = tsetf{2};
    spkData{ik} = tsetx.xset;
    fprintf(num2str(ik) + "-");

end

%% Load all spk condition specific

for ik = 1:Nfiles

    tsetx = load(spkpath + spkfiles{ik});
    tsetf = strsplit(spkfiles{ik}, "_");
    spkArea{ik} = tsetf{2};
    chidx = tsetx.xset.chids;

    for icond = 1:12

        sspkData{ik, icond} = tsetx.xset.xs{icond}(:, chidx, :);
        fprintf(num2str(icond) + "-");

    end

    fprintf(">" + num2str(ik) + "\n");

end

%% Load lfp

for ik = 1:Nfiles

    tsetx = load(lfppath + lfpfiles{ik});
    lfpData{ik} = tsetx.xset;
    fprintf(num2str(ik) + "-");

end

%% Labels

areaList = ["V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"];
aaab = 1;
axab = 2;
aaxb = 3;
aaax = 4;
bbba = 5;
bxba = 6;
bbxa = 7;
bbbx = 8;
rrrr = 9;
rxrr = 10;
rrxr = 11;
rrrx = 12;

ncntn = 0;
areaIDs = zeros(1, 1);
fileIDs = zeros(1, 1);

for ik = 1:Nfiles
        
    disp(spkArea{ik});
    areaIDtemp = find(strcmpi(spkArea{ik}, areaList));
    ncnt = size(sspkData{ik, 1}, 2);
    areaIDs(ncntn+1:ncntn+ncnt) = areaIDtemp;
    fileIDs(ncntn+1:ncntn+ncnt) = ik;
    disp(size(sspkData{ik, 1}));
    ncntn = ncntn + ncnt;

end

disp(ncntn);

%% Grand matrix concatenation (PEV) Stim Identity (A?B)

gkernel = ones(1, 1, 30); 

tOi1 = 501:4500;
tOi2 = 501:4500;
tN = length(tOi1);
icond1 = aaab;
icond2 = bbba;
nTrials = 200;

gmatrix1 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, icond1};
    tcnt = size(tempSig1, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = tempSig1(iTrials, :, tOi1);

    tempSig2 = sspkData{ik, icond2};
    ncnt = size(tempSig2, 2);
    tcnt = size(tempSig2, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = tempSig2(iTrials, :, tOi2);

    data = zeros(nTrials*2, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1, [1, 2], 1);
    gmatrix1(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(expv.*(p < 0.01));
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Grand matrix concatenation (PEV) Omission Identity (X|A?X|B)

tOi1 = 501:1500;
tOi2 = 1531:2530;
tOi3 = 2561:3560;
tOi4 = 3591:4590;

tN = length(tOi1);
% icond1 = axab;
% icond2 = aaxb;
nTrials = 100;

gmatrix2 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, axab};
    tcnt = size(tempSig1, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = tempSig1(iTrials, :, tOi2);

    tempSig2 = sspkData{ik, icond2};
    ncnt = size(tempSig2, 2);
    tcnt = size(tempSig2, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = tempSig2(iTrials, :, 3501:4500);

    data = zeros(nTrials*2, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrix2(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(expv);
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Grand matrix concatenation (iFR)

tN = 500;
icond1 = rrrr;
icond2 = rrrx;
nTrials = 100;

gmatrix3 = zeros(1, tN);
gmatrix4 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, icond1};
    % tcnt = size(tempSig1, 1);
    % iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = sum(tempSig1(:, :, 3501:4000), 1);
   
    tempSig2 = sspkData{ik, icond2};
    ncnt = size(tempSig2, 2);
    % tcnt = size(tempSig2, 1);
    % iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = sum(tempSig2(:, :, 3501:4000), 1);

    disp(ik);
    disp(length(unique(tempSig1)));
    disp(length(unique(tempSig2)));
    gmatrix3(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(tempSig1);
    gmatrix4(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(tempSig2);
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Spike scatter / Information / Information
% PEV(x^) | PEV(s>) (2Dscatter, color on area,)

% sPEVtime = 1:4500;
% xPEVtime = 1:1000;

sPEVs = zeros(1, size(gmatrix1, 1));
xPEVs = zeros(1, size(gmatrix2, 1));
sAFR = zeros(1, size(gmatrix3, 1));
xAFR = zeros(1, size(gmatrix4, 1));

for iN = 1:size(gmatrix1, 1)
    sPEVs(iN) = 100*(median(smooth(gmatrix1(iN, :), 10)));
    xPEVs(iN) = 100*(median(smooth(gmatrix2(iN, :), 10)));
    sAFR(iN) = median(gmatrix3(iN, :));
    xAFR(iN) = median(gmatrix4(iN, :));
end

colmap = ones(neuronCnt, 1);

figure;

subplot(2, 1, 1);
for iA = 1:11

    if iA < 5
        color_t = [1 0 0];
    elseif iA < 9
        color_t = [0 1 0];
    elseif iA < 10
        color_t = [0 0 1];
    else
        color_t = [0 0 0];
    end

    scatter(sPEVs(areaIDs == iA), xPEVs(areaIDs == iA), colmap(areaIDs == iA)*5, colmap(areaIDs == iA)*color_t, "filled", DisplayName=areaList(iA));
    % xlim([-1 1]);
    hold("on");
end
xlabel("Stim-PEV");ylabel("OXM-PEV");
legend;

subplot(2, 1, 2);
for iA = 1:11

    if iA < 5
        color_t = [1 0 0];
    elseif iA < 9
        color_t = [0 1 0];
    elseif iA < 10
        color_t = [0 0 1];
    else
        color_t = [0 0 0];
    end

    scatter(sAFR(areaIDs == iA), xAFR(areaIDs == iA), colmap(areaIDs == iA)*5, colmap(areaIDs == iA)*color_t, "filled", DisplayName=areaList(iA));
    xlim([0 20]);
    ylim([0 20]);
    line(1:20, 1:20);
    hold("on");
end
xlabel("Stim-AFR");ylabel("OXM-AFR");
legend;

%%



%%

% {cond}[trial,neuron,time]
% PCA<
% Kmeans<

% xrefs = 
% xrefx = 

isx04 = sspkData{3, 2};
isx08 = sspkData{3, 10};

imx1 = squeeze(mean(isx04, 1));
imagesc(imx1);

%%
temp_sig1 = isx04(:, 50, :);
temp_sig2 = isx08(:, 50, :);
temp_sig1 = smooth(squeeze(mean(temp_sig1, 1)), 50);
temp_sig2 = smooth(squeeze(mean(temp_sig2, 1)), 50);
plot(temp_sig1 - mean(temp_sig1));
hold("on");
plot(temp_sig2 - mean(temp_sig2));

%%
mspkx = isx09;
mspkxo = mspkx;
mspkxo = smoothdata2(mspkxo, "gaussian", {1, 20});
mspkx = smoothdata2(mspkx, "gaussian", {1, 40});

kg = 12;
[i1, c1, sm1, d1] = kmeans(mspkx, kg, "Distance", "sqeuclidean");

figure;
sgtitle("All neurons>XA?XB>PEV>Kmeans");

for ik = 1:kg
    subplot(kg, 1, ik);
    plot(100*mean(mspkxo(i1 == ik, :)));
    xline(500, "Color", "#DD0000");
    xline(1530, "Color", "#DD0000");
    xline(2560, "Color", "#DD0000");
    xline(3590, "Color", "#DD0000");
    cntx = sum(i1 == ik);
    ylabel("PEV%");
    title("Kgroup" + num2str(ik) + ", nNeuron = " + num2str(cntx));
end

%%

% --- Script to Plot a Random Image in Multiple Tabs ---

% Clear workspace, command window, and close all figures
clear;
clc;
close all;

% --- Configuration ---
imageSize = [100, 100]; % Define the dimensions [height, width] of the image
numberOfTabs = 4;      % Define how many tabs you want in the figure

% --- Generate Random Image ---
% Create a 100x100 matrix with random values between 0 and 1
randomImage = rand(imageSize(1), imageSize(2));

% --- Create Multi-Tab Figure ---
% Create a new figure window (using uifigure for modern UI components)
fig = uifigure('Name', 'Multi-Tab Random Image Display', 'Position', [100 100 600 500]);

% Create a tab group that fills most of the figure window
tabGroup = uitabgroup(fig, 'Position', [20 20 fig.Position(3)-40 fig.Position(4)-40]);

% --- Populate Tabs with the Image ---
% Loop through the desired number of tabs
for i = 1:numberOfTabs
    % Create a new tab within the tab group
    tab = uitab(tabGroup, 'Title', ['Tab ' num2str(i)]);

    % Create an axes object within the current tab
    % Using uiaxes is recommended with uifigure and uitab
    ax = uiaxes(tab);

    % Display the random image in the current axes
    % imagesc scales the data to use the full colormap, which is often
    % useful for visualizing data.
    imagesc(ax, randomImage);

    % Set the colormap to grayscale (common for this type of data)
    colormap(ax, gray);

    % Ensure the image has the correct aspect ratio (pixels are square)
    axis(ax, 'image');

    % Turn off the axis ticks for a cleaner image view (optional)
    ax.XTick = [];
    ax.YTick = [];

    % Add a title to the plot within the tab
    title(ax, ['Random 100x100 Image - Plot ' num2str(i)]);

    % Add a colorbar to show the mapping of data values to colors
    colorbar(ax);
end

% --- End of Script ---

%%

%% Bench

% FEF : 44-45-53 (+> 19-70 )

isx09 = sspkData{1, 9};
isx10 = sspkData{1, 10};

im1 = squeeze(mean(isx09(:, :, 1:4000), 1));
im2 = squeeze(mean(isx10(:, :, 1:4000), 1));

for ik = 1:size(im1, 1)

    % im1(ik, :) = (im1(ik, :) - mean(im1(ik, :))) / std(im1(ik, :));
    % im2(ik, :) = (im2(ik, :) - mean(im2(ik, :))) / std(im2(ik, :));
    im1(ik, :) = smooth(im1(ik, :), 50);
    im2(ik, :) = smooth(im2(ik, :), 50);

end

figure;
subplot(2, 2, 1);
imagesc(im1);
% clim([-2 2]);
subplot(2, 2, 2);
imagesc(im2);
% clim([-2 2]);
subplot(2, 2, 3);
imagesc(im1 - im2);
% clim([-2 2]);
subplot(2, 2, 4);
imagesc(im1 ./ im2);
% clim([-5 5]);

%%
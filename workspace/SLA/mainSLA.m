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
spkData = cell(Nfiles, 1);
sspkData = cell(Nfiles, 12);

%% LFP unifier

areax = "FEF";
lfppath = "lfpSet\";
lfpfiles = {dir(lfppath).name};
lfpfiles = lfpfiles(contains(lfpfiles, areax));
Nfiles = length(lfpfiles);
lfpData = cell(Nfiles, 1);

%% Load spk

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

for ik = 1:Nfiles
        
    disp(spkArea{ik});
    disp(size(sspkData{ik, 1}));

end

%% Grand matrix concatenation 

tN = size(sspkData{1, 1}, 3);
icond1 = aaab;
icond2 = bbba;

gmatrix1 = zeros(250, 1, tN);
gmatrix2 = zeros(250, 1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig = sspkData{ik, icond1};
    ncnt = size(tempSig, 2);
    tcnt = size(tempSig, 1);
    iTrials = mod(randperm(250), tcnt) + 1;
    gmatrix1(:, neuronCnt+1:neuronCnt+ncnt, :) = tempSig(iTrials, :, :);

    tempSig = sspkData{ik, icond2};
    ncnt = size(tempSig, 2);
    tcnt = size(tempSig, 1);
    iTrials = mod(randperm(250), tcnt) + 1;
    gmatrix2(:, neuronCnt+1:neuronCnt+ncnt, :) = tempSig(iTrials, :, :);

    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Grand PEVs



%% Spike scatter / Information / Information
% data = zeros(N*2, n1, 1500);
% data(1:N, :, :) = x1;
% data(N+1:2*N, :, :) = x2;
% % data(2*N+1:end, :, :) = x3;
% 
% % data = jSmooth(data, 100);
% % data(:, :, 1:500) = data(:, :, randperm(500, 500));
% 
% groupIDs = [ones(1, N), ones(1, N)*2];%, ones(1, N)*3];
% [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
% expvard{k}(ntotal+1:ntotal+n1, :) = squeeze(expv);
% {cond}[trial,neuron,time]
% PCA<
% Kmeans<
% PEV(x^) | PEV(s>) (2Dscatter, color on area)
% xrefs = 
% xrefx = 

ia = 1;
icond = 8;
isx09 = sspkData{ia, icond};
disp(spkArea{ia});

mspkx = squeeze(sum(isx09, 1));
mspkxo = mspkx;
mspkxo = smoothdata2(mspkxo, "gaussian", {1, 20});
mspkx = smoothdata2(mspkx, "gaussian", {1, 40});

kg = 7;
[i1, c1, sm1, d1] = kmeans(mspkx, kg, "Distance", "sqeuclidean");

figure;
sgtitle(spkArea{ia} + ">cond>" + num2str(icond));

for ik = 1:kg
    subplot(kg, 1, ik);
    plot(mean(mspkxo(i1 == ik, :)));
    xline(500, "Color", "#DD0000");
    xline(1530, "Color", "#DD0000");
    xline(2560, "Color", "#DD0000");
    xline(3590, "Color", "#DD0000");
    cntx = sum(i1 == ik);
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
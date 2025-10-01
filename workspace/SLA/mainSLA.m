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
spkData = cell(Nfiles, 1);
sspkData1 = cell(Nfiles, 1);
sspkData2 = cell(Nfiles, 1);

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
    spkData{ik} = tsetx.xset;
    fprintf(num2str(ik));

end

%% Load spk specific condition

scond1 = 9;
scond2 = 10;

for ik = 1:Nfiles

    tsetx = load(spkpath + spkfiles{ik});
    sspkData1{ik} = tsetx.xset.xs{scond1};
    sspkData2{ik} = tsetx.xset.xs{scond2};
    fprintf(num2str(ik));

end

%% Load lfp

for ik = 1:Nfiles

    tsetx = load(lfppath + lfpfiles{ik});
    lfpData{ik} = tsetx.xset;
    fprintf(num2str(ik));

end

%% Bench

% FEF : 44-45-53 (+> 19-70 )

isx09 = spkData{1}.xs{9};
isx10 = spkData{1}.xs{10};

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

%% Spike scatter

% xrefs = 
% xrefx = 
isx09 = spkData{1}.xs{7};

mspkx = squeeze(sum(isx09, 1));
mspkx = smoothdata2(mspkx, "gaussian", {1, 20});

kg = 3;
[i1, c1, sm1, d1] = kmeans(mspkx, kg, "Distance", "sqeuclidean");

figure;
for ik = 1:kg
    subplot(kg, 1, ik);
    plot(mean(mspkx(i1 == ik, :)));
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
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

%% TFR unifier

areax = "V1";
tfrpath = "tfrSet\";
tfrfiles = {dir(tfrpath).name};
tfrfiles = tfrfiles(contains(tfrfiles, areax));
Nfiles = length(tfrfiles);
tfrData = cell(Nfiles, 1);

%%

for ik = 1:Nfiles

    tfrx = load(tfrpath + tfrfiles{ik});
    tfrData{ik} = tfrx.tset.pgx;
    fprintf(num2str(ik));

end

%% Bench

jTFRLaminarPlotter(tfrData{4}, tsinfo, "V1", 50, 40, 20);

%% Deep vs Sup
%% Alignment
%% Stemplots
%% Functions

function jTFRLaminarPlotter(tfrdata, tset, areaname, chn, chd, cl4)

    im1 = tfrdata{3} + tfrdata{7} + tfrdata{11};
    
    tbaselinex = tset.tbands{1}(5:end-5);
    
    for ik = 1:size(im1, 2)
    
        for jk = 1:size(im1, 1)
    
            im1(jk, ik, :) = im1(jk, ik, :) / mean(im1(jk, ik, tbaselinex), "all");
    
        end
    
    end
    
    fmapx = tset.fmap;
    locx = (linspace(1, chn, chn) - cl4) * chd;
    
    figure;
    
    tctx1 = tset.tbands{3}(end-20:end);
    imx1 = 10*log(squeeze(mean(im1(:, :, tctx1), 3)));
    imx1 = smoothdata2(imx1, "movmedian", 10);
    
    subplot(2, 2, 1);
    imagesc(imx1, "XData", fmapx, "YData", locx);
    yline(0);
    xlabel("Freq.");
    ylabel("Dist. from L4 in um");
    title(areaname + " (baseline before omission)");
    clim([-15 15]);
    % set(gca, "YDir", "normal");
    cb = colorbar();
    ylabel(cb, "Power vs. baseline (dB)");
    
    tctx2 = tset.tbands{4}(1:30);
    imx2 = 10*log(squeeze(mean(im1(:, :, tctx2), 3)));
    imx2 = smoothdata2(imx2, "movmedian", 10);
    
    subplot(2, 2, 2);
    imagesc(imx2, "XData", fmapx, "YData", locx);
    yline(0);
    xlabel("Freq.");
    ylabel("Dist. from L4 in um");
    title(areaname + " (omission)");
    clim([-15 15]);
    % set(gca, "YDir", "normal");
    cb = colorbar();
    ylabel(cb, "Power vs. baseline (dB)");
    
    subplot(2, 1, 2);
    imagesc(imx2 - imx1, "XData", fmapx, "YData", locx);
    yline(0);
    xlabel("Freq.");
    ylabel("Dist. from L4 in um");
    title(areaname + " (omission - pre-omission-base)");
    clim([-10 10]);
    % set(gca, "YDir", "normal");
    cb = colorbar();
    ylabel(cb, "Power vs. baseline (dB)");

end

%%

figure;
% subplot(2, 1, 1);
stem(fmapx, mean(imx3(1:41, :), 1), "DisplayName", "Deep(?)");
% title("Sup");
% subplot(2, 1, 2);
hold("on");
stem(fmapx, mean(imx3(41:end, :), 1), "DisplayName", "Sup(?)");
legend();
title("S(after omission) vs. S(first stim fx)");

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
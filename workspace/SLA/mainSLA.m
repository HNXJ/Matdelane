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

% V1/V2/V3d/V3a

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

jTFRLaminarPlotter(tfrData{2}, tsinfo, "V1", 52, 40, 16, 1);

%% Deep vs Sup concat by area

ls = [16, 15, 22, 20, 20]; % V1
lflip = [0, 0, 1, 0, 1]; % V1

% l4s = [21, 22, 20, 22]; % V4
% lflip = [0, 0, 0, 1];

% l4s = [15, 20, 4, 13, 8]; % MT
% lflip = [1, 0, 1, 1, 1];

% l4s = [50, 10, 80, 40, 30]; % PFC
% lflip = [1, 1, 0, 1, 0];

tdxsup = cell(12, 1);
tdxdeep = cell(12, 1);

for ik = 1:Nfiles

    if lflip(ik)

        for jk = 1:12
    
            tdxdeep{jk} = [tdxdeep{jk}; tfrData{ik}{jk}(1:l4s(ik), :, :)];
            tdxsup{jk} = [tdxsup{jk}; tfrData{ik}{jk}(l4s(ik)+1:end, :, :)];
    
        end

    else

        for jk = 1:12
    
            tdxsup{jk} = [tdxsup{jk}; tfrData{ik}{jk}(1:l4s(ik), :, :)];
            tdxdeep{jk} = [tdxdeep{jk}; tfrData{ik}{jk}(l4s(ik)+1:end, :, :)];
    
        end

    end

end


%% Stemplots

[a1, a2, a3] = jTFRLaminarPlotter(tdxsup, tsinfo, "Vx-sup", 50, 40, 20);
[b1, b2, b3] = jTFRLaminarPlotter(tdxdeep, tsinfo, "Vx-deep", 50, 40, 20);

%%
% 1. Finish population TFRx results
% 
% 2. Statistical info for TFRx
% 
% 3. Position (S3|X2 , S4|X3) vs (S3|f , S4|f) control
% 
% 4. Spectral shifts (TFR : stim-base vs stim-base-ox)
% 
% 4.1 Characterize center frequency (FOOOF/iP)
%%

figure;
subplot(2, 2, 1); % Spectrolam
subplot(2, 2, 2); % Stem allCh Sox vs Sfx
subplot(2, 2, 3); % Stem deepCh/supCh (+-2SEM|chN)
subplot(2, 2, 4); % Barplot deepCh-supCh X AlphaBeta<8-30>|Gamma<32-150>

%%

fmapx = tsinfo.fmap;
deepx = mean(a3, 1);
supx = mean(b3, 1);

figure;
subplot(2, 1, 1);
stem(fmapx, deepx, "DisplayName", "Deep");

hold("on");
stem(fmapx, supx, "DisplayName", "Sup");
legend();
title("Stim after omission vs stim after fixation power change");

subplot(2, 1, 2);
stem(fmapx, deepx-supx);
title("S(after omission) vs. S(first stim fx) deep + vs sup -");

%% Functions

function [imx1x, imx2x, imx3x] = jTFRLaminarPlotter(tfrdata, tset, areaname, chn, chd, cl4, vflipx)

    im1 = tfrdata{2} + tfrdata{6} + tfrdata{10};
    % im1 = tfrdata{3} + tfrdata{7} + tfrdata{11};
    
    tbaselinex = tset.tbands{2}(end-20:end);

    if ~exist("vflipx", "var")

        vflipx = 0;

    end

    if vflipx

        for ik = 1:size(im1, 2)
        
            im1(:, ik, :) = im1(:, ik, :) / mean(im1(:, ik, tbaselinex), "all");
        
        end

    else

        for ik = 1:size(im1, 2)
    
            for jk = 1:size(im1, 1)
    
                im1(jk, ik, :) = im1(jk, ik, :) / mean(im1(jk, ik, tbaselinex), "all");
    
            end
    
        end

    end
    
    fmapx = tset.fmap;
    locx = (linspace(1, chn, chn) - cl4) * chd;
    
    figure;
    
    tctx1 = tset.tbands{4};
    imx1x = squeeze(mean(im1(:, :, tctx1), 3));
    imx1 = 10*log(imx1x);
    % imx1 = smoothdata2(imx1, "movmedian", 10);
    
    subplot(2, 2, 1);
    imagesc(imx1, "XData", fmapx, "YData", locx);
    yline(0);
    xticks(0:20:200)
    xlabel("Freq.");
    ylabel("Dist. um");
    title(areaname + " (S after omission)");
    clim([-10 10]);
    xlim([0 150]);
    % set(gca, "YDir", "normal");
    cb = colorbar();
    ylabel(cb, "Power vs. baseline (dB)");
    
    tctx2 = tset.tbands{2};
    imx2x = squeeze(mean(im1(:, :, tctx2), 3));
    imx2 = 10*log(imx2x);
    imx2 = smoothdata2(imx2, "movmedian", 10);
    
    subplot(2, 2, 2);
    imagesc(imx2, "XData", fmapx, "YData", locx);
    yline(0);
    xticks(0:20:200)
    xlabel("Freq.");
    ylabel("Dist. um");
    title(areaname + " (S after fixation)");
    clim([-10 10]);
    xlim([0 150]);
    % set(gca, "YDir", "normal");
    cb = colorbar();
    ylabel(cb, "Power vs. baseline (dB)");
    
    imx3x = 100*(imx1x - imx2x);
    subplot(2, 1, 2);
    imagesc(imx3x, "XData", fmapx, "YData", locx);
    yline(0);
    xticks(0:20:200)
    xlabel("Freq.");
    ylabel("Dist. from L4 in um");
    title(areaname + " (S-x vs S-o)");
    % clim([-50 50]);
    xlim([0 150]);
    % set(gca, "YDir", "normal");
    cb = colorbar();
    ylabel(cb, "Change (%)");

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
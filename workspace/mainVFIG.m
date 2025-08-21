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

areax = "PFC"; % V1-V2-V3d-V3a-V4-MT-MST-TEO-FST-FEF-PFC
tfrpath = "tfrSet\";
tfrfiles = {dir(tfrpath).name};
tfrfiles = tfrfiles(contains(tfrfiles, areax));
Nfiles = length(tfrfiles);
tfrData = cell(Nfiles, 1);

%

for ik = 1:Nfiles

    tfrx = load(tfrpath + tfrfiles{ik});
    tfrData{ik} = tfrx.tset.pgx;
    fprintf(num2str(ik));

end

%

tfrsx = cell(5, 12);
foidx = {[11:28], [32:48], [52:100], [140:320], [321:800]};

for ik = 1:Nfiles

    for jk = 1:5

        for kk = 1:12
        
            tfrsx_temp = tfrData{ik}{kk}(:, foidx{jk}, :);

            if ik == 1
            
                tfrsx{jk, kk} = tfrsx_temp;

            else

                tfrsx{jk, kk} = cat(1, tfrsx{jk, kk}, tfrsx_temp);

            end

        end

    end

    disp(ik);

end

save("tfrSetBand\" + areax, "tfrsx", "-v7.3");

%% TFR coherence (CorrXspearman)

tfrpath = "tfrSetBand\";
tfrfiles = {dir(tfrpath).name};
tfrfiles = tfrfiles(contains(tfrfiles, ".mat"));
Nfiles = length(tfrfiles);
tfrData = cell(Nfiles, 1);

%%

for ik = 1:Nfiles

    tfrx = load(tfrpath + tfrfiles{ik});
    tfrData{ik} = tfrx.tfrsx;
    fprintf(num2str(ik));

end

%%

% toix = tsinfo.tbands{2}(1:45);
% toix = tsinfo.tbands{3}(16:60);
% toix = 1:299;
cond = 10;

for fband = 1:5

    gcorr{fband} = zeros(1, length(toix));
    lcnt = 0;

    areacnt = zeros(1, Nfiles);
    
    for ik = 1:Nfiles
    
        ncnt = size(tfrData{ik}{fband, cond}, 1);
        % nbands = size(tfrData{ik}{fband, cond}, 2);
    
        temp_gmat = squeeze(mean(tfrData{ik}{fband, cond}(:, :, toix), 2));
        
        for jk = 1:ncnt

            temp_gmat(jk, :) = temp_gmat(jk, :) / mean(temp_gmat(jk, end-10:end));

        end

        sub_temp_gmat = sum(corr(temp_gmat', "Type", "Spearman"));
        [~, indxsorted] = sort(sub_temp_gmat, "descend");
        temp_gmat = temp_gmat(indxsorted, :);
        % temp_gmat = smoothdata2(temp_gmat, "movmean", {51 15});

        gcorr{fband}(lcnt+1:lcnt+ncnt, :) = temp_gmat;
    
        lcnt = lcnt + ncnt;
        areacnt(ik) = lcnt;

    end

end

% gcorr_stim = gcorr;
% gcorr_omxn = gcorr;

%%


fbandlabels = {"Theta[2.5-7Hz]", "Alpha[8-12Hz]", "Beta[13-30Hz]", "L-Gamma[32-80Hz]", "H-Gamma[80-200Hz]"};
figure;

for fbandx = 1:4

    imx1 = corr(gcorr{fbandx}', "Type", "Spearman");
    % imx1 = smoothdata2(imx1, "movmean", 51);
    imx1(abs(imx1) < 0.5);

    subplot(2, 2, fbandx);

    imagesc(imx1);
    
    for ik = 1:Nfiles
    
        xline(areacnt(ik), "Color", [1 0 0], "LineWidth", 1);
        yline(areacnt(ik), "Color", [1 0 0], "LineWidth", 1);
    
    end
    
    title("Spectral power corr" + fbandlabels{fbandx})
    colormap("parula");

end

sgtitle("Omission time xcorr"); % TODO make |r| < .25 0

%%

fbandlabels = {"Theta[2.5-7Hz]", "Alpha[8-12Hz]", "Beta[13-30Hz]", "L-Gamma[32-80Hz]", "H-Gamma[80-200Hz]"};
areanames = {"V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"};

figure("Position", [0 0 1200 2200]);

for fbandx = 1:4

    imx1 = corr(gcorr_stim{fbandx}', "Type", "Spearman");

    gcorr_a = zeros(11, 11);

    for ik = 1:Nfiles

        if ik == 1

            lcntx = 0;

        else

            lcntx = areacnt(ik-1);

        end
        
        for jk = 1:ik

            if jk == 1
    
                rcntx = 0;
    
            else
    
                rcntx = areacnt(jk-1);
    
            end

            gcorr_a(ik, jk) = mean(imx1(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");

        end

    end

    % imx1(abs(imx1) < 0.5);

    subplot(2, 2, fbandx);

    gcorr_a(abs(gcorr_a) < 0.01) = 0;
    gcorr_a = gcorr_a*100;

    imagesc(gcorr_a);

    for ik = 1:Nfiles

        for jk = 1:ik
            
            text(jk-.3, ik, num2str(gcorr_a(ik, jk), 3))
            
        end

    end
    clim([-100 100]);
    
    title("Spectral power corr" + fbandlabels{fbandx});
    xticks(1:11);
    xticklabels(areanames);
    yticks(1:11);
    yticklabels(areanames);
    colormap("parula");

end


sgtitle("Sx time Corr "); % TODO make |r| < .25 0
fname = "allareatfr_corr" + "stim";
print(gcf,'-vector','-dsvg', fname +".svg");

%%

fbandlabels = {"Theta[2.5-7Hz]", "Alpha[8-12Hz]", "Beta[13-30Hz]", "L-Gamma[32-80Hz]", "H-Gamma[80-200Hz]"};
areanames = {"V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"};

figure("Position", [0 0 1200 2200]);

for fbandx = 1:4

    imx1 = corr(gcorr_stim{fbandx}', "Type", "Spearman");
    imx2 = corr(gcorr_omxn{fbandx}', "Type", "Spearman");

    gcorr_a = zeros(11, 11);
    gcorr_b = zeros(11, 11);

    for ik = 1:Nfiles

        if ik == 1

            lcntx = 0;

        else

            lcntx = areacnt(ik-1);

        end
        
        for jk = 1:ik

            if jk == 1
    
                rcntx = 0;
    
            else
    
                rcntx = areacnt(jk-1);
    
            end

            gcorr_a(ik, jk) = mean(imx1(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");
            gcorr_b(ik, jk) = mean(imx2(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");

        end

    end

    % imx1(abs(imx1) < 0.5);

    subplot(2, 2, fbandx);

    gcorr_c = 100*(gcorr_b.^2 - gcorr_a.^2);
    gcorr_c(isnan(gcorr_c)) = 0;
    gcorr_c(abs(gcorr_c) < 1) = 0;

    imagesc(gcorr_c);

    for ik = 1:Nfiles

        for jk = 1:ik

            if gcorr_c(ik, jk)
                
                text(jk-.3, ik, num2str(gcorr_c(ik, jk), 3));

            else

                text(jk-.3, ik, "n.s");

            end

        end

    end
    clim([-20 20]);
    
    title("Spectral power corr ratio stim to omission" + fbandlabels{fbandx});
    xticks(1:11);
    xticklabels(areanames);
    yticks(1:11);
    yticklabels(areanames);
    colormap("parula");

end

sgtitle("R2 change (omxn - stim)"); % TODO make |r| < .25 0
fname = "allareatfr_corr" + "_r2change";
print(gcf,'-vector','-dsvg', fname +".svg");

%% Barplot

fbandlabels = {"Theta[2.5-7Hz]", "Alpha[8-12Hz]", "Beta[13-30Hz]", "L-Gamma[32-80Hz]", "H-Gamma[80-200Hz]"};
areanames = {"V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"};

% figure("Position", [0 0 1200 2200]);

gcorr_g = zeros(11, 11);

for fbandx = 5:5

    imx1 = corr(gcorr_stim{fbandx}', "Type", "Spearman");
    imx2 = corr(gcorr_omxn{fbandx}', "Type", "Spearman");

    gcorr_a = zeros(11, 11);
    gcorr_b = zeros(11, 11);

    for ik = 1:Nfiles

        if ik == 1

            lcntx = 0;

        else

            lcntx = areacnt(ik-1);

        end
        
        for jk = 1:Nfiles

            if jk == 1
    
                rcntx = 0;
    
            else
    
                rcntx = areacnt(jk-1);
    
            end

            gcorr_a(ik, jk) = mean(imx1(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");
            gcorr_b(ik, jk) = mean(imx2(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");

        end

    end

    % imx1(abs(imx1) < 0.5);

    % subplot(2, 2, fbandx);

    gcorr_c = 100*(gcorr_a.^2 + gcorr_b.^2);
    gcorr_c(isnan(gcorr_c)) = 0;
    gcorr_c(abs(gcorr_c) < 1) = 0;

    gcorr_g = gcorr_g + gcorr_c;
    % 
    % bar(mean(gcorr_a, 1));
    % 
    % % for ik = 1:Nfiles
    % % 
    % %     for jk = 1:ik
    % % 
    % %         text(jk-.3, ik, num2str(gcorr_a(ik, jk), 3))
    % % 
    % %     end
    % % 
    % % end
    % ylim([0 30]);
    % 
    % title("Spectral power corr" + fbandlabels{fbandx});
    % xticks(1:11);
    % xticklabels(areanames);
    % % yticks(1:11);
    % % yticklabels(areanames);
    % colormap("parula");

end
 
figure;

barrsq = mean(abs(gcorr_g));

bar(barrsq);
xticks(1:11);
xticklabels(areanames);

sgtitle("Sx time Corr "); % TODO make |r| < .25 0
fname = "allareatfr_corr" + "stimbar";
print(gcf,'-vector','-dsvg', fname +".svg");
%%
a = [1024 384 128 384 384 512 640 512 512 1152 1280]

%% Bench

sigtemp = tfrData{4};
disp(size(sigtemp{1}));
nChlx = size(sigtemp{1}, 1);
dSpacing = 40;
l4Ch = 12;
jTFRLaminarPlotter(sigtemp, tsinfo, areax, nChlx, dSpacing, l4Ch, 1);

%%

immx = tfrData{3}{2};
immx = squeeze(mean(immx, 1));
figure;imagesc(log(immx), "XData", linspace(0, 200, 801));

%% Deep vs Sup concat by area

l4s = [17, 15, 22, 10, 12]; % V1
lflip = [0, 0, 0, 0, 0]; % V1

% l4s = [22, 12, 30, 28]; % V2
% lflip = [0, 0, 1, 1]; % V2

% l4s = [25, 30, 22, 40, 16, 30]; % V3d
% lflip = [1, 0, 1, 0, 0, 0]; % V3d

% l4s = [18, 18, 30]; % V3a
% lflip = [0, 1, 1]; % V3a

% l4s = [30, 30, 28, 32]; % V4
% lflip = [1, 0, 0, 0];

% l4s = [20, 30, 4, 36, 2]; % MT
% lflip = [0, 0, 1, 0, 0];

% l4s = [49, 25, 20, 24, 18]; % MST
% lflip = [0, 0, 0, 0, 1];

% l4s = [40, 40]; % TEO
% lflip = [1, 1];

% l4s = 20; % FST
% lflip = 1;

% l4s = [48, 25, 27]; % FEF
% lflip = [0, 1, 0];

% l4s = [60, 50, 10, 60, 30]; % PFC
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
deepxse = std(a3, 1, 1) / sqrt(size(tdxsup{1}, 1));
supxse = std(b3, 1, 1) / sqrt(size(tdxdeep{1}, 1));
supx = mean(b3, 1);

cl1 = [0.3 0.3 0.9];
cl2 = [0.9 0.3 0.3];

figure;
subplot(2, 1, 1);
plot(fmapx, deepx, "DisplayName", "Deep");
patch([fmapx; fmapx(end:-1:1)], [deepx + 2*deepxse, deepx(end:-1:1) - 2*deepxse(end:-1:1)], cl1, "EdgeColor", "none", "FaceColor", cl1, "FaceAlpha", 0.5, "HandleVisibility", "off");
yline(0, "--", "HandleVisibility", "off");
hold("on");
plot(fmapx, supx, "DisplayName", "Sup");
xlabel("Frequency (Hz)");
ylabel("Percentage change");
patch([fmapx; fmapx(end:-1:1)], [supx + 2*supxse, supx(end:-1:1) - 2*supxse(end:-1:1)], cl2, "EdgeColor", "none", "FaceColor", cl2, "FaceAlpha", 0.5, "HandleVisibility", "off");
legend();
title("Stim after omission vs stim after fixation power change");

diffx = deepx - supx;
subplot(2, 1, 2);
plot(fmapx, diffx);
xlabel("Frequency (Hz)");
ylabel("Percentage change");
yline(0, "--");
patch([fmapx; fmapx(end:-1:1)], [diffx + deepxse + supxse, diffx(end:-1:1) - deepxse(end:-1:1) - supxse(end:-1:1)], cl1, "EdgeColor", "none", "FaceColor", cl1, "FaceAlpha", 0.5, "HandleVisibility", "off");
title("S(after omission) vs. S(first stim fx) deep + vs sup -");

sgtitle("Population; +-2SE; (nCh ~= 500) | " + areax);

%% Functions

function [imx1x, imx2x, imx3x] = jTFRLaminarPlotter(tfrdata, tset, areaname, chn, chd, cl4, vflipx)

    im1 = tfrdata{2} + tfrdata{6} + tfrdata{10};
    % im1 = tfrdata{3} + tfrdata{7} + tfrdata{11};
    
    tbaselinex = tset.tbands{1};

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
    imx1 = smoothdata2(imx1, "movmedian", 10);
    
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
    imx3x = smoothdata2(imx3x, "movmedian", 10);
    subplot(2, 1, 2);
    imagesc(imx3x, "XData", fmapx, "YData", locx);
    yline(0);
    xticks(0:20:200)
    xlabel("Freq.");
    ylabel("Dist. from L4 in um");
    title(areaname + " (S-x vs S-o)");
    clim([-50 50]);
    xlim([0 150]);
    % set(gca, "YDir", "normal");
    cb = colorbar();
    ylabel(cb, "Change (%)");

end

%%


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
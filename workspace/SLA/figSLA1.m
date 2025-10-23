%% Setup

color_t = zeros(12, 3);
color_t(1, :) = [.75 0 0];
color_t(2, :) = [1 0 0];
color_t(3, :) = [.8 .4 0];
color_t(4, :) = [.75 .75 0];
color_t(5, :) = [.4 .8 0];
color_t(6, :) = [0 1 0];
color_t(7, :) = [0 .8 .4];
color_t(8, :) = [.1 .3 .8];
color_t(9, :) = [0 0 1];
color_t(10, :) = [.3 0 .7];
color_t(11, :) = [.45 0 .55];
color_t(12, :) = [.15 .15 .75];

%% Scatter vectors

sPEV = zeros(1, size(gmatrix1, 1));
xPEV = zeros(1, size(gmatrix2, 1));
xOPEV = zeros(1, size(gmatrix6, 1));
sAFR = zeros(1, size(gmatrix3, 1));
xAFR = zeros(1, size(gmatrix4, 1));
bAFR = zeros(1, size(gmatrix5, 1));
sVFR = zeros(1, size(gmatrix3, 1));
xVFR = zeros(1, size(gmatrix4, 1));
bVFR = zeros(1, size(gmatrix5, 1));

for iN = 1:size(gmatrix1, 1)
    sPEV(iN) = 100*(max(smooth(gmatrix1(iN, :), 100)));
    xPEV(iN) = 100*(max(smooth(gmatrix2(iN, :), 100)));
    xOPEV(iN) = 100*(max(smooth(gmatrix6(iN, :), 100)));

    svPEV(iN) = 100*(std(smooth(gmatrix1(iN, :), 50)));
    xvPEV(iN) = 100*(std(smooth(gmatrix2(iN, :), 50)));
    tW = size(gmatrix3, 2) / 1000;
    sAFR(iN) = sum(gmatrix3(iN, :)) / tW;
    xAFR(iN) = sum(gmatrix4(iN, :)) / tW;
    sVFR(iN) = var(gmatrix3(iN, :)) / tW;
    xVFR(iN) = var(gmatrix4(iN, :)) / tW;
    tW = size(gmatrix5, 2) / 1000;
    bAFR(iN) = sum(gmatrix5(iN, :)) / tW;
    bVFR(iN) = var(gmatrix5(iN, :)) / tW;
end

colmap = ones(neuronCnt, 1);

%% Fig.1 A:
%  aFR(Oxm) | aFR(Stim)

figure("Position", [0 0 1500 1500]);

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sAFR(areaIDs == iA);
    ye1 = xAFR(areaIDs == iA);
    be1 = bAFR(areaIDs == iA);
    idxs = be1 > 0.1;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), 1, color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .6, 'std');
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.1, "HandleVisibility", "off", "EdgeColor", "none");
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.0, "HandleVisibility", "off", "EdgeColor", color_t(iA, :), "LineStyle", ":", "LineWidth", 2);

    xetext = max(xe);
    yetext = max(ye);
    text(xetext, yetext, [areaList{iA}, ''], "Color", color_t(iA, :), "FontWeight", "bold");

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .6, 'std');
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.0, "HandleVisibility", "off", "EdgeColor", color_t(iA, :), "LineStyle", ":", "LineWidth", 2);
   
    line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    
end

set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'linear');
xlim([0.3 30]);
ylim([0.3 30]);
% xlim([0.1 20]);
% ylim([0.1 20]);
grid("on");
% zlim([0.2 5000]);
set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

title("1.a: Omission|Stimulus avg. FR");
xlabel("Stim-AFR");ylabel("Oxm-AFR");
legend;

% fname = "1A-Scatter-Ox-St";
% print(gcf,'-vector','-dsvg', fname +".svg");

%% Fig.1 B:
%  Relative to baseline scatter aFR(Oxm/b) | aFR(Stim/b)

figure("Position", [0 0 1500 1500]);

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    be1 = bAFR(areaIDs == iA);
    xe1 = sAFR(areaIDs == iA)./be1;
    ye1 = xAFR(areaIDs == iA)./be1;
    idxs = be1 > 0.1;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), 1, color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    idxs = be1 > 1.0 | xe1 > 1.0;
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, 1, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.05, "HandleVisibility", "off", "EdgeColor", "none");
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.0, "HandleVisibility", "off", "EdgeColor", color_t(iA, :), "LineStyle", ":", "LineWidth", 1);

    xetext = mean(xe);
    yetext = mean(ye);
    text(xetext, yetext, [areaList{iA}, ''], "Color", color_t(iA, :), "FontWeight", "bold");

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.1, "HandleVisibility", "off", "EdgeColor", "none");
   
    line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");

    
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([0.1 100]);
ylim([0.1 100]);

set(gca, 'XTick', [0.1 1 2 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 2 3 10 30 100]);
set(gca, 'YTick', [0.1 1 2 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 2 3 10 30 100]);

title("1.b: Omission/baseline|Stimulus/baseline avg. FR ratio");
xlabel("Stim-AFR/baseline");ylabel("Oxm-AFR/baseline");
legend;

% fname = "1B-Scatter-Ox-St-BaselineRatio";
% print(gcf,'-vector','-dsvg', fname +".svg");

%% Fig.1 C:
%  Baseline scatter aFR(Oxm) | aFR(baseline)

figure;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = bAFR(areaIDs == iA);
    ye1 = xAFR(areaIDs == iA);
    ze1 = sAFR(areaIDs == iA);
    idxs = xe1 > .1;
    idxs2 = ye1 ./ xe1 > 2 & xe1 > 1;
    pointsizes = 3*ones(size(generalIDs));
    pointsizes(idxs2) = 40;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), pointsizes(idxs), color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
        ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, 1, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.05, "HandleVisibility", "off", "EdgeColor", "none");

    xetext = mean(xe);
    yetext = mean(ye);
    text(xetext, yetext, [areaList{iA}, ''], "Color", color_t(iA, :), "FontWeight", "bold");

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.1, "HandleVisibility", "off", "EdgeColor", "none");
   
    line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([0.1 30]);
ylim([0.1 30]);
set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

title("1.c: Omission|Baseline avg. FR");
xlabel("Base-AFR");ylabel("Oxm-AFR");
legend;

%% Fig.1 D:
%  Kmeans clustered aFR(Oxm) | aFR(Stim)

figure;

kmeans_gn = 7;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sAFR(areaIDs == iA); % gm3
    ye1 = xAFR(areaIDs == iA); % gm4
    be1 = sAFR(areaIDs == iA); % gm5
    ze1 = bVFR(areaIDs == iA);
    we1 = xVFR(areaIDs == iA);

    pointsizes = 3*ones(size(generalIDs));

    kGx = [xe1; ye1; ze1; we1]';
    [idxs, c1, sm1, d1] = kmeans(kGx, kmeans_gn, "Distance", "sqeuclidean");
    
    for iB = 1:kmeans_gn

        idxs2 = find(idxs == iB);
        scatter3(xe1(idxs2), ye1(idxs2), idxs2, pointsizes(idxs2), color_t(iA, :), "filled", DisplayName=areaList(iA) + " G." + num2str(iB));
        view(0, 90);
        hold("on");
        
        if length(idxs2) > 3

                ze = mean(generalIDs);
                [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, 1, 'std');
                patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");
            
                xetext = mean(xe);
                yetext = mean(ye);
                text(xetext, yetext, [areaList{iA}, ''], "Color", color_t(iA, :), "FontWeight", "bold");

                [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, .5, 'std');
                patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.4, "HandleVisibility", "off", "EdgeColor", "none");
               
                line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
                xlim([0.2 100]);
                ylim([0.2 100]);
        end

    end

    line(0:100, 0:100, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    xlim([0.3 30]);
    ylim([0.3 30]);
    
end

xlim([0.3 100]);
ylim([0.3 100]);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

title("1.d: Omission|Stimulus avg. FR (Kmeans groups = " + num2str(kmeans_gn) + ")");
xlabel("Stim-AFR");ylabel("Oxm-AFR");
legend;

%% Fig.1 E:
%  Kmeans clustered aFR(Oxm) | aFR(Baseline)

figure;

kmeans_gn = 7;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = bAFR(areaIDs == iA); % gm3
    ye1 = xAFR(areaIDs == iA); % gm4
    be1 = sAFR(areaIDs == iA); % gm5
    ze1 = bVFR(areaIDs == iA);
    we1 = xVFR(areaIDs == iA);

    pointsizes = 3*ones(size(generalIDs));

    kGx = [xe1; ye1; ze1; we1]';
    [idxs, c1, sm1, d1] = kmeans(kGx, kmeans_gn, "Distance", "cityblock");
    
    for iB = 1:kmeans_gn

        idxs2 = find(idxs == iB);
        scatter3(xe1(idxs2), ye1(idxs2), idxs2, pointsizes(idxs2), color_t(iA, :), "filled", DisplayName=areaList(iA) + " G." + num2str(iB));
        view(0, 90);
        hold("on");
        
        if length(idxs2) > 3

            ze = mean(generalIDs);
            [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, 1, 'std');
            patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");
        
            xetext = mean(xe);
            yetext = mean(ye);
            text(xetext, yetext, [areaList{iA}, ''], "Color", color_t(iA, :), "FontWeight", "bold");

            [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, .5, 'std');
            patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.4, "HandleVisibility", "off", "EdgeColor", "none");
           
            line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
            xlim([0.2 100]);
            ylim([0.2 100]);

        end

    end

    line(0:100, 0:100, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");

    
end

xlim([0.3 100]);
ylim([0.3 100]);
set(gca, 'XScale', 'log', 'YScale', 'log');

title("1.e: Omission|Baseline avg. FR (Kmeans groups = " + num2str(kmeans_gn) + ")");
xlabel("Base-AFR");ylabel("Oxm-AFR");
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

legend;

%% Scatter variation(Oxm) | variation (Stim)

%% Fig.2 A:
%  Baseline scatter PEV(oxm) | PEV(stim)

figure;

xrks = cell(2, 11);

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sPEV(areaIDs == iA);
    ye1 = xPEV(areaIDs == iA);
    ze1 = bAFR(areaIDs == iA);
    idxs = xe1 > .5 | ye1 > .5;
    idxs2 = ye1 ./ xe1 > 0 | ye1 > 1;
    pointsizes = 1*ones(size(generalIDs));
    pointsizes(idxs2) = 5;

    xrks{1, iA} = xe1(idxs);
    xrks{2, iA} = ye1(idxs);

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), pointsizes(idxs), color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .6, 'std');
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.1, "HandleVisibility", "off", "EdgeColor", "none");
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.0, "HandleVisibility", "off", "EdgeColor", color_t(iA, :), "LineStyle", ":", "LineWidth", 2);

    xetext = max(xe);
    yetext = max(ye);
    text(xetext, yetext, [areaList{iA}, ''], "Color", color_t(iA, :), "FontWeight", "bold");

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .6, 'std');
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.0, "HandleVisibility", "off", "EdgeColor", color_t(iA, :), "LineStyle", ":", "LineWidth", 2);
   
    % ze = mean(generalIDs);
    % [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, 1, 'std');
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.1, "HandleVisibility", "off", "EdgeColor", "none");
    % 
    % xetext = mean(xe);
    % yetext = mean(ye);
    % text(xetext, yetext, [areaList{iA}, ''], "Color", color_t(iA, :), "FontWeight", "bold");
    % 
    % [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, .5, 'std');
    % patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");
   
    line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([0.1 100]);
ylim([0.1 100]);

set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

title("2.a: Omission|Stim PEV");
xlabel("Stim-PEV");ylabel("Oxm-PEV");
legend;

%% Fig.2 b:
%  Baseline scatter PEV(oxm) | AFR(oxm/base)

figure;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = xAFR(areaIDs == iA) ./ bAFR(areaIDs == iA);
    ye1 = xPEV(areaIDs == iA);
    ze1 = bAFR(areaIDs == iA);
    idxs = ze1 > 0.1;
    idxs2 = ye1 > 2;
    pointsizes = 1*ones(size(generalIDs));
    pointsizes(idxs2) = 2;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), pointsizes(idxs), color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, 1, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");

    re = - 0.5 * ((-1)^ceil(iA/2))* xe + ((-1)^iA)*ye;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.4, "HandleVisibility", "off", "EdgeColor", "none");
   
    line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    xlim([0.1 30]);
    ylim([0.1 30]);
    
end

set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

title("2.b: Omission PEV|Omission/Baseline avg. FR");
xlabel("Oxm-AFR/Base-AFR");ylabel("Oxm-PEV");
legend;

%% Fig.2 c:
%  PEV(oxm_pos) | AFR(oxm/base)

figure;
title("2.c: Omission position PEV|Omission/Baseline avg. FR");
xlabel("Oxm-AFR/Base-AFR");ylabel("OxmPos-PEV");
legend;

%% Fig.2 d:
%  PEV(oxm_pos) | AFR(oxm/base)

figure;
title("2.d: Omission position PEV|Omission identity PEV");
xlabel("Oxm-PEV");ylabel("OxmPos-PEV");
legend;

%% Fig.2 E:
%  Kmeans clustered PEV(Oxm) | PEV(Stim)

figure;

kmeans_gn = 3;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sPEV(areaIDs == iA); % gm3
    ye1 = xPEV(areaIDs == iA); % gm4
    be1 = sAFR(areaIDs == iA); % gm5
    ze1 = bVFR(areaIDs == iA);
    we1 = xVFR(areaIDs == iA);

    pointsizes = 3*ones(size(generalIDs));

    kGx = [xe1; ye1; ze1; we1]';
    [idxs, c1, sm1, d1] = kmeans(kGx, kmeans_gn, "Distance", "cityblock");
    
    for iB = 1:kmeans_gn

        idxs2 = find(idxs == iB);
        scatter3(xe1(idxs2), ye1(idxs2), idxs2, pointsizes(idxs2), color_t(iA, :), "filled", DisplayName=areaList(iA) + " G." + num2str(iB));
        view(0, 90);
        hold("on");
        
        if length(idxs2) > 3

            ze = mean(generalIDs);
            [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, 1, 'std');
            patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");
        
            re = - 0.5 * ((-1)^ceil(iA/2))* xe + ((-1)^iA)*ye;
            idtextz = find(re == max(re));
            text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));
        
            [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, .5, 'std');
            patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.4, "HandleVisibility", "off", "EdgeColor", "none");
           
            line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
            xlim([0.1 100]);
            ylim([0.1 100]);

        end

    end

    line(0:100, 0:100, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");

    
end

xlim([0.1 100]);
ylim([0.1 100]);
set(gca, 'XScale', 'log', 'YScale', 'log');

title("2.e: Omission|Stim PEV (Kmeans groups = " + num2str(kmeans_gn) + ")");
xlabel("Stim-PEV");ylabel("Oxm-PEV");

set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

legend;

%% Fig.2 F:
%  Kmeans clustered PEV(Oxm) | PEV(Oxm-Pos)

figure;

kmeans_gn = 3;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = xOPEV(areaIDs == iA); % gm3
    ye1 = xPEV(areaIDs == iA); % gm4
    be1 = sAFR(areaIDs == iA); % gm5
    ze1 = bVFR(areaIDs == iA);
    we1 = xvPEV(areaIDs == iA);

    pointsizes = 3*ones(size(generalIDs));

    kGx = [xe1; ye1; ze1; we1]';
    [idxs, c1, sm1, d1] = kmeans(kGx, kmeans_gn, "Distance", "cityblock");
    
    for iB = 1:kmeans_gn

        idxs2 = find(idxs == iB);
        scatter3(xe1(idxs2), ye1(idxs2), idxs2, pointsizes(idxs2), color_t(iA, :), "filled", DisplayName=areaList(iA) + " G." + num2str(iB));
        view(0, 90);
        hold("on");
        
        if length(idxs2) > 3

            ze = mean(generalIDs);
            [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, 1, 'std');
            patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", "none");
        
            re = - 0.5 * ((-1)^ceil(iA/2))* xe + ((-1)^iA)*ye;
            idtextz = find(re == max(re));
            text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));
        
            [xe, ye] = fitConfidenceEllipse(xe1(idxs2), ye1(idxs2), 1000, .5, 'std');
            patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.4, "HandleVisibility", "off", "EdgeColor", "none");
           
            line([0.1 100], [0.1 100], "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
            xlim([0.1 100]);
            ylim([0.1 100]);

        end

    end

    line(0:100, 0:100, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");

    
end

xlim([0.1 100]);
ylim([0.1 100]);
set(gca, 'XScale', 'log', 'YScale', 'log');

title("2.e: Omission|Omission Pos-PEV (Kmeans groups = " + num2str(kmeans_gn) + ")");
xlabel("Oxm-Pos-PEV");ylabel("Oxm-PEV");

set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTick', [0.1 1 3 10 30 100]);
set(gca, 'XTickLabel', [0.1 1 3 10 30 100]);
set(gca, 'YTick', [0.1 1 3 10 30 100]);
set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);

legend;

%% Functions

function [X, Y] = generateEllipse(x, y, a, b, theta, N)
    N = max(3, round(N));

    t = linspace(0, 2*pi, N);

    x_prime = a * cos(t);
    y_prime = b * sin(t);
    
    cos_t = cos(theta);
    sin_t = sin(theta);

    x_rotated = x_prime .* cos_t - y_prime .* sin_t;
    y_rotated = x_prime .* sin_t + y_prime .* cos_t;

    X = x_rotated + x;
    Y = y_rotated + y;

end

function [X_ellipse, Y_ellipse] = fitConfidenceEllipse(X_data, Y_data, N, CI, mode, rect)
    
    if ~exist('mode', 'var')

        mode = 'sem';

    end

    if ~exist('rect', 'var')

        rect = 1;

    end

    % wx = sum(X_data(:));wy = sum(Y_data(:));
    % dataMatrix = [X_data(:).^2/wx, Y_data(:).^2/wy];
    dataMatrix = [X_data(:), Y_data(:)];
    
    mu = mean(dataMatrix);
    center_x = mu(1);
    center_y = mu(2);
    
    C = cov(dataMatrix);
    
    [V, D] = eig(C);
    
    eigenvalues = diag(D);
    
    [lambda_sorted, idx] = sort(eigenvalues, 'descend');
    
    lambda_a = lambda_sorted(1);
    lambda_b = lambda_sorted(2);
    
    V_major = V(:, idx(1));
    theta = atan2(V_major(2), V_major(1));

    if strcmpi(mode, 'chi2')

        DEGREES_OF_FREEDOM = 2;
        chi2_critical_value = chi2inv(CI, DEGREES_OF_FREEDOM);
        k = sqrt(chi2_critical_value); 

    elseif strcmpi(mode, 'std')

        k = CI; 

    elseif strcmpi(mode, 'sem')
        
        k = CI / sqrt(length(X_data));

    else

        k = 1.0;

    end

    radius_a = k * sqrt(lambda_a);
    radius_b = k * sqrt(lambda_b);
    
    [X_ellipse, Y_ellipse] = generateEllipse(center_x, center_y, radius_a, radius_b, theta, N);

    if rect

        X_ellipse(X_ellipse < 0) = 1e-6;
        Y_ellipse(Y_ellipse < 0) = 1e-6;
        
    end

end

%% Other options

%  Sorted Raster (image, 4099neuron in 4500ms, sorted)
%  Scatter PEV(Oxm) | PEV(Stim) ; Time-Concat version (to resolve offset / presence issue)
%  Scatter PEV(Oxm-A?B) | PEV(Stim)
%  Scatter PEV(Oxm-2?3?4) | PEV(Stim)
%  Scatter PEV(Stim after omission) | PEV(Stim)

%% Ranksums

disp("Wilcoxon rank sum test.");
disp("Oxm vs Stim (PEV vectors)");
figure;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sPEV(areaIDs == iA);
    ye1 = xPEV(areaIDs == iA);
    ze1 = bAFR(areaIDs == iA);
    idxs = ze1 > .1;

    xrks{1, iA} = xe1(idxs);
    xrks{2, iA} = ye1(idxs);

    [p, h, stt] = ranksum(xrks{1, iA}, xrks{2, iA}, 'tail', 'right');
    disp(areaList(iA) + " ; Mean(Stm,Oxm) = (" + num2str(mean(xrks{1, iA})) + "," + num2str(mean(xrks{2, iA})) + ")");
    scatter(mean(xrks{1, iA}), mean(xrks{2, iA}), 40, color_t(iA, :), "filled");hold("on")
    if p < 0.01
        text(1+mean(xrks{1, iA}), mean(xrks{2, iA}), "Stim->" + num2str(p, 1));
    else
        text(1+mean(xrks{1, iA}), mean(xrks{2, iA}), "n.s");
    end
    text(mean(xrks{1, iA}), mean(xrks{2, iA}), areaList{iA});
    disp(" p = " + num2str(p) + " / H0-reject(Skewed-right?) : " + num2str(h));

    [p, h, stt] = ranksum(xrks{1, iA}, xrks{2, iA}, 'tail', 'left');
    if p < 0.01
        text(mean(xrks{1, iA}), 0.1+mean(xrks{2, iA}), "Oxm-^" + num2str(p, 1));
    else
        text(mean(xrks{1, iA}), 0.1+mean(xrks{2, iA}), "n.s");
    end
    disp(" p = " + num2str(p) + " / H0-reject(Skewed-left?) : " + num2str(h));

end

xlim([0 20]);
ylim([1 3]);
%%

for ikA = 1:5

    ar1 = ikA*2 - 1;
    ar2 = ikA*2;
    figure;
    
    subplot(2, 1, 1);
    histogram(xrks{1, ar1}, 80, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Stim");    
    hold("on");
    xline(mean(xrks{1, ar1}), "Color", "blue");

    histogram(xrks{2, ar1}, 80, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Oxm");
    xline(mean(xrks{2, ar1}), "Color", "red");
    title("PEV-PDF" + areaList(ar1));
    
    subplot(2, 1, 2);
    histogram(xrks{1, ar2}, 80, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Stim");
    hold("on");
    xline(mean(xrks{1, ar2}), "Color", "blue");
    histogram(xrks{2, ar2}, 80, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Oxm");
    xline(mean(xrks{2, ar2}), "Color", "red");
    title("PEV-PDF" + areaList(ar2));
    legend;

end

%%

ar1 = 1;
ar2 = 11;
figure;

subplot(2, 1, 1);
histogram(xrks{1, ar1}, 100, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Stim");
hold("on");
histogram(xrks{2, ar1}, 100, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Oxm");
title("PEV-PDF" + areaList(ar1));

subplot(2, 1, 2);
histogram(xrks{1, ar2}, 100, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Stim");
hold("on");
histogram(xrks{2, ar2}, 100, "BinLimits", [0 20], "Normalization", "probability", "DisplayName", "Oxm");
title("PEV-PDF" + areaList(ar2));
legend;

    %%
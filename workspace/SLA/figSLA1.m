%% Setup

color_t = zeros(11, 3);
color_t(1, :) = [1 .2 0];
color_t(2, :) = [.9 .4 0];
color_t(3, :) = [.8 .6 0];
color_t(4, :) = [.6 .8 0];
color_t(5, :) = [.4 .7 0];
color_t(6, :) = [.2 .8 0];
color_t(7, :) = [0 .6 .6];
color_t(8, :) = [0 .4 .8];
color_t(9, :) = [.3 .2 .7];
color_t(10, :) = [.6 .1 .8];
color_t(11, :) = [.5 0 1];

%% Scatter vectors

sPEV = zeros(1, size(gmatrix1, 1));
xPEV = zeros(1, size(gmatrix2, 1));
sAFR = zeros(1, size(gmatrix3, 1));
xAFR = zeros(1, size(gmatrix4, 1));
bAFR = zeros(1, size(gmatrix5, 1));
sVFR = zeros(1, size(gmatrix3, 1));
xVFR = zeros(1, size(gmatrix4, 1));
bVFR = zeros(1, size(gmatrix5, 1));

for iN = 1:size(gmatrix1, 1)
    sPEV(iN) = 100*(max(smooth(gmatrix1(iN, :), 50)));
    xPEV(iN) = 100*(max(smooth(gmatrix2(iN, :), 50)));
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

figure;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sAFR(areaIDs == iA);
    ye1 = xAFR(areaIDs == iA);
    be1 = bAFR(areaIDs == iA);
    idxs = be1 > 1.0;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), 1, color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .2, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.6, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.3, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
   
    line(0:20, 0:20, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    xlim([0 20]);
    ylim([0 20]);
    
end

title("1.a: Omission|Stimulus avg. FR");
xlabel("Stim-AFR");ylabel("Oxm-AFR");
legend;

%% Fig.1 B:
%  Relative to baseline scatter aFR(Oxm/b) | aFR(Stim/b)

figure;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    be1 = bAFR(areaIDs == iA);
    xe1 = sAFR(areaIDs == iA)./be1;
    ye1 = xAFR(areaIDs == iA)./be1;
    idxs = be1 > 1.0;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), 1, color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .2, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.6, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.3, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
   
    line(0:20, 0:20, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    xlim([0 4]);
    ylim([0 4]);
    
end

title("1.b: Omission/baseline|Stimulus/baseline avg. FR ratio");
xlabel("Stim-AFR/baseline");ylabel("Oxm-AFR/baseline");
legend;

%% Fig.1 C:
%  Baseline scatter aFR(Oxm) | aFR(baseline)

figure;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = bAFR(areaIDs == iA);
    ye1 = xAFR(areaIDs == iA);
    ze1 = sAFR(areaIDs == iA);
    idxs = xe1 > 1 | ye1 > 1;
    idxs2 = ye1 ./ xe1 > 2;
    pointsizes = 3*ones(size(generalIDs));
    pointsizes(idxs2) = 40;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), pointsizes(idxs), color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .2, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.6, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.3, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
   
    line(1:20, 1:20, "color", [0 0 0], "HandleVisibility", "off");
    xlim([0 20]);
    ylim([0 20]);
    
end

title("1.c: Omission|Baseline avg. FR");
xlabel("Base-AFR");ylabel("Oxm-AFR");
legend;

%% Fig.1 D:
%  Kmeans clustered aFR(Oxm) | aFR(Stim)

figure;

kmeans_gn = 4;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sAFR(areaIDs == iA);
    ye1 = xAFR(areaIDs == iA);
    be1 = bAFR(areaIDs == iA);

    % [idxs, c1, sm1, d1] = kmeans(xe1, kmeans_gn, "Distance", "sqeuclidean");
    % TODO for each group / elliptic fit
    scatter(xe1(idxs), ye1(idxs), 1, color_t(iA, :), "filled", DisplayName=areaList(iA));
    hold("on");
    
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .2, 'std');
    patch(xe, ye, color_t(iA, :), "FaceAlpha", 0.3, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off", "EdgeColor", [1 1 1]);

    line(0:20, 0:20, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    xlim([0 20]);
    ylim([0 20]);
    
end

title("1.d: Omission|Stimulus avg. FR (Kmeans groups = " + num2str(kmeans_gn) + ")");
xlabel("Stim-AFR");ylabel("Oxm-AFR");
legend;

%% Scatter variation(Oxm) | variation (Stim)

%% Fig.2 A:
%  Baseline scatter PEV(oxm) | PEV(stim)

figure;

for iA = 1:11

    generalIDs = find(areaIDs == iA);
    xe1 = sPEV(areaIDs == iA);
    ye1 = xPEV(areaIDs == iA);
    ze1 = bAFR(areaIDs == iA);
    idxs = xe1 > 1 | ye1 > 1;
    idxs2 = ye1 ./ xe1 > 0 & ye1 > 2.19;
    pointsizes = 1*ones(size(generalIDs));
    pointsizes(idxs2) = 10;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), pointsizes(idxs), color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .2, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.6, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.3, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
   
    line(0:20, 0:20, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    xlim([0 20]);
    ylim([0 5]);
    
end

title("2.a: Omission|Baseline avg. FR");
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
    pointsizes(idxs2) = 10;

    scatter3(xe1(idxs), ye1(idxs), generalIDs(idxs), pointsizes(idxs), color_t(iA, :), "filled", DisplayName=areaList(iA));
    view(0, 90)
    hold("on");
    
    ze = mean(generalIDs);
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .2, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.6, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ''], "Color", color_t(iA, :));

    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, .5, 'std');
    patch(xe, ye, ze*ones(size(xe)), color_t(iA, :), "FaceAlpha", 0.3, "HandleVisibility", "off", "EdgeColor", [1 1 1]);
   
    line(0:20, 0:20, "color", [0 0 0], "HandleVisibility", "off", "LineStyle", "--");
    xlim([0 5]);
    ylim([0 5]);
    
end

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

function [X_ellipse, Y_ellipse] = fitConfidenceEllipse(X_data, Y_data, N, CI, mode)
    
    if ~exist('mode', 'var')

        mode = 'sem';

    end

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

end

%% Other options

%  Sorted Raster (image, 4099neuron in 4500ms, sorted)
%  Scatter PEV(Oxm) | PEV(Stim) ; Time-Concat version (to resolve offset / presence issue)
%  Scatter PEV(Oxm-A?B) | PEV(Stim)
%  Scatter PEV(Oxm-2?3?4) | PEV(Stim)
%  Scatter PEV(Stim after omission) | PEV(Stim)

%%
%%

color_t = zeros(11, 3);
color_t(1, :) = [1 .2 0];
color_t(2, :) = [.9 .4 0];
color_t(3, :) = [.8 .6 0];
color_t(4, :) = [.6 .8 0];
color_t(5, :) = [.4 1 0];
color_t(6, :) = [.2 .8 0];
color_t(7, :) = [0 .6 .6];
color_t(8, :) = [0 .4 .8];
color_t(9, :) = [.3 .2 .7];
color_t(10, :) = [.6 0 .8];
color_t(11, :) = [.9 0 1];

%% Spike scatter / Information / Information

sPEVs = zeros(1, size(gmatrix1, 1));
xPEVs = zeros(1, size(gmatrix2, 1));
sAFR = zeros(1, size(gmatrix3, 1));
xAFR = zeros(1, size(gmatrix4, 1));
bAFR = zeros(1, size(gmatrix5, 1));

for iN = 1:size(gmatrix1, 1)
    sPEVs(iN) = 100*(max(smooth(gmatrix1(iN, :), 10)));
    xPEVs(iN) = 100*(max(smooth(gmatrix2(iN, :), 10)));
    tW = size(gmatrix3, 2) / 1000;
    sAFR(iN) = sum(gmatrix3(iN, :)) / tW;
    xAFR(iN) = sum(gmatrix4(iN, :)) / tW;
    tW = size(gmatrix5, 2) / 1000;
    bAFR(iN) = sum(gmatrix5(iN, :)) / tW;
end

colmap = ones(neuronCnt, 1);

%% Fig.1 panels:
%  Scatter aFR(Oxm) | aFR(Stim)

figure;

for iA = 9:11

    xe1 = sAFR(areaIDs == iA);
    ye1 = xAFR(areaIDs == iA);
    be1 = bAFR(areaIDs == iA);
    idxs = be1 > 2.0;

    scatter(xe1(idxs), ye1(idxs), 1, color_t(iA, :), "filled", DisplayName=areaList(iA));
    hold("on");
    
    [xe, ye] = fitConfidenceEllipse(xe1(idxs), ye1(idxs), 1000, 0.5);
    patch(xe, ye, color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off");
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), [areaList{iA}, ' - S+'], "Color", color_t(iA, :));

    % [xe, ye] = fitConfidenceEllipse(xe1, ye1, 1000, 0.5);
    % patch(xe, ye, color_t(iA, :), "FaceAlpha", 0.2, "HandleVisibility", "off");
    % re = xe.^2 + ye.^2;
    % idtextz = find(re == max(re));
    % text(xe(idtextz), ye(idtextz), [areaList{iA}, ' - all'], "Color", color_t(iA, :));

    xlim([0 20]);
    ylim([0 20]);
    
end

xlabel("Stim-AFR");ylabel("Oxm-AFR");
legend;
%%
figure;

for iA = 1:11

    scatter(sAFR(areaIDs == iA), xAFR(areaIDs == iA), colmap(areaIDs == iA)*1, color_t(iA, :), "filled", DisplayName=areaList(iA));
    hold("on");
    center_p = [mean(sAFR(areaIDs == iA)), mean(xAFR(areaIDs == iA))];
    dF1 = sum(areaIDs == iA);
    a = std(sPEVs(areaIDs == iA));%/ sqrt(dF1);
    b = std(xPEVs(areaIDs == iA));%/ sqrt(dF1);
    theta = atan(mean(xAFR(areaIDs == iA) > 0.1) / mean(sAFR(areaIDs == iA) > 0.1));
    [xe, ye] = generateEllipse(center_p(1), center_p(2), a, b, theta, 1000);
    plot(xe, ye, "LineWidth", 1, "Color", color_t(iA, :), "HandleVisibility", "off");
    re = xe.^2 + ye.^2;
    idtextz = find(re == max(re));
    text(xe(idtextz), ye(idtextz), areaList{iA}, "Color", color_t(iA, :));
    xlim([0 20]);
    ylim([0 20]);
    line(1:20, 1:20, "HandleVisibility", "off");
    hold("on");

end

xlabel("Stim-AFR");ylabel("OXM-AFR");
legend;

%%  Scatter variation(Oxm) | variation (Stim)
%  Sorted Raster (image, 4099neuron in 4500ms, sorted)


%%

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

function [X_ellipse, Y_ellipse] = fitConfidenceEllipse(X_data, Y_data, N, k)
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
    
    radius_a = k * sqrt(lambda_a);
    radius_b = k * sqrt(lambda_b);
    
    [X_ellipse, Y_ellipse] = generateEllipse(center_x, center_y, radius_a, radius_b, theta, N);

end

% function [X_ellipse, Y_ellipse] = fitConfidenceEllipse(X_data, Y_data, N)
% 
%     dataMatrix = [X_data(:), Y_data(:)];
% 
%     mu = mean(dataMatrix);
%     center_x = mu(1);
%     center_y = mu(2);
% 
%     C = cov(dataMatrix);
% 
%     [V, D] = eig(C);
% 
%     eigenvalues = diag(D);
% 
%     [lambda_sorted, idx] = sort(eigenvalues, 'descend');
% 
%     lambda_a = lambda_sorted(1);
%     lambda_b = lambda_sorted(2);
% 
%     V_major = V(:, idx(1));
%     theta = atan2(V_major(2), V_major(1)); 
% 
%     CONFIDENCE_LEVEL = 0.20;
%     DEGREES_OF_FREEDOM = 2;
% 
%     chi2_critical_value = chi2inv(CONFIDENCE_LEVEL, DEGREES_OF_FREEDOM);
%     k = sqrt(chi2_critical_value);
% 
%     radius_a = k * sqrt(lambda_a);
%     radius_b = k * sqrt(lambda_b);
% 
%     [X_ellipse, Y_ellipse] = generateEllipse(center_x, center_y, radius_a, radius_b, theta, N);
% 
% end

%%

% Fig.2 panels:
%  Scatter PEV(Oxm) | PEV(Stim) ; Time-Concat version (to resolve offset / presence issue)
%  Scatter PEV(Oxm-A?B) | PEV(Stim)
%  Scatter PEV(Oxm-2?3?4) | PEV(Stim)
%  Scatter PEV(Stim after omission) | PEV(Stim)


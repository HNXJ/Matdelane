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
    sPEV(iN) = 100*(max(smoothdata(gmatrix1(iN, :), "gaussian", 100)));
    xPEV(iN) = 100*(max(smoothdata(gmatrix2(iN, :), "gaussian", 100)));
    xOPEV(iN) = 100*(max(smoothdata(gmatrix6(iN, :), "gaussian", 100)));

    svPEV(iN) = 100*(std(smoothdata(gmatrix1(iN, :), "gaussian", 100)));
    xvPEV(iN) = 100*(std(smoothdata(gmatrix2(iN, :), "gaussian", 100)));
    tW = size(gmatrix3, 2) / 1000;
    sAFR(iN) = sum(gmatrix3(iN, :)) / tW;
    sVFR(iN) = var(gmatrix3(iN, :)) / tW;
    tW = size(gmatrix4, 2) / 1000;
    xAFR(iN) = sum(gmatrix4(iN, :)) / tW;
    xVFR(iN) = var(gmatrix4(iN, :)) / tW;
    tW = size(gmatrix5, 2) / 1000;
    bAFR(iN) = sum(gmatrix5(iN, :)) / tW;
    bVFR(iN) = var(gmatrix5(iN, :)) / tW;
end

colmap = ones(neuronCnt, 1);

%% Stim prime neurons

g1_sprime = find(sAFR./bAFR > 2.5);
g2_invsprime = find(xAFR./sAFR > 1.5);
g3_inv_sprime = find(sAFR < bAFR*3 & bAFR > 5.0);

%% Oxm prime neurons

g1_oxprime = find(xPEV > 0);

%% Group PEVs in time (N) with iFRs

nIDgroup = g1_sprime;
% nIDgroup = g2_invsprime;
% nIDgroup = g1_oxprime;

% nIDgroup = g2_inv_sprime;

kW = 500;
kX = 2;
tN = 1000; % length(temp_sig1);
timevec = linspace(-500, 4250, 4750);

figure("Position", [0 0 1500 1500]);

areaNcountx = zeros(1, 11);

for iA = 1:11

    cntxtemp = sum(areaIDs(nIDgroup) == iA);
    areaNcountx(iA) = cntxtemp;

end
subplot(2, 1, 1);
bar(areaNcountx);
xticks(1:11);
xticklabels(areaList(1:11));

condOi = [1, 2, 9];
ncondOi = length(condOi);

subplot(2, 1, 2);

for ik = 1:ncondOi
    
    icond = condOi(ik);
    temp_sig1 = zeros(size(timevec));
    temp_sig2 = zeros(size(timevec));

    for nIDx = nIDgroup

        temp_sigxc = find(squeeze(sspkDataCleanTrials{fileIDs(nIDx), icond}(:, infileIDs(nIDx))) == 1);
        temp_sigx = squeeze(sspkData{fileIDs(nIDx), icond}(temp_sigxc, infileIDs(nIDx), :));
        temp_sigx = tN*smoothdata2(temp_sigx, "gaussian", {1, kW});
        temp_sigx(isnan(temp_sigx)) = 0;
        
        if size(temp_sigx, 1) > 1
            temp_sig1 = temp_sig1 + squeeze(mean(temp_sigx, 1))/length(nIDgroup);
            temp_sig2 = temp_sig2 + squeeze(std(temp_sigx, 1) / sqrt(2*size(temp_sigx, 1)))/(2*length(nIDgroup));
        end

    end

    % temp_sig1 = detrend(temp_sig1);

    xpe1 = timevec;
    ype1 = temp_sig1 + kX*temp_sig2;
    xpe2 = timevec(end:-1:1);
    ype2 = temp_sig1(end:-1:1) - kX*temp_sig2(end:-1:1);
    % plot(timevec, temp_sig1, "LineWidth", 2, "Color", color_t(icond, :));
    patch([xpe1, xpe2], [ype1, ype2], color_t(icond, :), "FaceAlpha", 0.5, "EdgeColor", "none", "DisplayName", condNames(icond) + "+-" + num2str(kX) + "SEM");
    hold("on");

end

xline(0, "HandleVisibility", "off");
xline(1030, "HandleVisibility", "off");
xline(2060, "HandleVisibility", "off");
xline(3090, "HandleVisibility", "off");

xlim([-750, 4500]);
legend();
xlabel("Time(ms)");ylabel("FR(Spk/s)");
timevec = linspace(-500, 4250, 4750);
dispnametemp1 = "PEV(AAXBvsBBXA)";
temp_sigx = gmatrixN2(nIDgroup, :);
for ik = 1:size(temp_sigx, 1)
    temp_sigx(ik, :) = temp_sigx(ik, :) - min(temp_sigx(ik, :));
end
temp_sig1 = squeeze(max(temp_sigx));
% temp_sig1 = temp_sig1 - min(temp_sig1);
yyaxis("right");

pevvec = 100*smoothdata(temp_sig1, "gaussian", kW);
pevvec = pevvec - min(pevvec);

% kW = 1;
plot(timevec, pevvec, "DisplayName", dispnametemp1, "LineWidth", 2, "color", [.5 0 .5]);
yline(0, "LineStyle", "--");

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ylabel("PEV(%)", "Color", [1 0 0]);
sgtitle("N = " + num2str(length(nIDgroup)));

%%
%% Spectral

[ps0, fs0, ts0] = pspectrum(smoothdata(temp_sig1(1:1000), "gaussian", 100), 1000, "power", "FrequencyLimits", [1 100]);

[ps1, fs1, ts1] = pspectrum(smoothdata(temp_sig1, "gaussian", 100), 1000, "power", "FrequencyLimits", [1 100]);

figure;
% imagsesc(log(ps1), "XData", ts1-0.5, "YData", fs1+0.1);
plot(fs1+0.1, 10*log10(ps1./ps0));
xlabel("Freq");
ylabel("change(dB)");
sgtitle("Neuron no." + num2str(nID) + " > " + areaList(areaIDs(nID)) + "/Power change from baseline");
% set(gca, "YDir", "normal");
set(gca, "XScale", "log");
set(gca, 'XTick', [1 3 8 12 20 40 100]);
set(gca, 'XTickLabel', [1 3 8 12 20 40 100]);
% set(gca, 'YTick', [0.1 1 3 10 30 100]);
% set(gca, 'YTickLabel', [0.1 1 3 10 30 100]);
%%
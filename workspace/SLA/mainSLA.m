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
condNames = ["aaab", "axab", "aaxb", "aaax", "bbba", "bxba", "bbxa", "bbbx", "rrrr", "rxrr", "rrxr", "rrrx"];

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
    infileIDs(ncntn+1:ncntn+ncnt) = 1:ncnt;
    disp(size(sspkData{ik, 1}));
    ncntn = ncntn + ncnt;

end

disp(ncntn);

%% Grand matrix concatenation (PEV) Omission Identity (A?B)

gkernel = ones(1, 1, 100); 

tOi1 = 1:4500;
tOi2 = 1:4500;
% tOi3 = 1:4500;
tN = length(tOi1);
icond1 = aaax;
icond2 = bbbx;
% icond3 = rrrx;
nTrials = 100;

gmatrixN1 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, icond1};
    ncnt = size(tempSig1, 2);
    tcnt = size(tempSig1, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = tempSig1(iTrials, :, tOi1);

    tempSig2 = sspkData{ik, icond2};
    tcnt = size(tempSig2, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = tempSig2(iTrials, :, tOi2);

    % tempSig3 = sspkData{ik, icond3};
    % tcnt = size(tempSig3, 1);
    % iTrials = mod(randperm(nTrials), tcnt) + 1;
    % tempSig3 = tempSig3(iTrials, :, tOi3);

    data = zeros(nTrials*2, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    % data(nTrials*2+1:3*nTrials, :, :) = tempSig3;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];%, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN1(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(expv);
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Grand matrix concatenation (PEV) Stim Identity (A?B)

gkernel = ones(1, 1, 100); 

tOi1 = 1:4500;
tOi2 = 1:4500;
tN = length(tOi1);
icond1 = aaab;
icond2 = bbba;
nTrials = 200;

gmatrixN2 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, icond1};
    ncnt = size(tempSig1, 2);
    tcnt = size(tempSig1, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = tempSig1(iTrials, :, tOi1);

    tempSig2 = sspkData{ik, icond2};
    tcnt = size(tempSig2, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = tempSig2(iTrials, :, tOi2);

    data = zeros(nTrials*2, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN2(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(expv);
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Grand matrix concatenation (PEV) Omission Identity (A?B)

gkernel = ones(1, 1, 100); 

tOi1 = 1:4500;
tOi2 = 1:4500;
% tOi3 = 1:4500;
tN = length(tOi1);
icond1 = axab;
icond2 = bxba;
% icond3 = rrrx;
nTrials = 100;

gmatrixN3 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, icond1};
    ncnt = size(tempSig1, 2);
    tcnt = size(tempSig1, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = tempSig1(iTrials, :, tOi1);

    tempSig2 = sspkData{ik, icond2};
    tcnt = size(tempSig2, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = tempSig2(iTrials, :, tOi2);

    % tempSig3 = sspkData{ik, icond3};
    % tcnt = size(tempSig3, 1);
    % iTrials = mod(randperm(nTrials), tcnt) + 1;
    % tempSig3 = tempSig3(iTrials, :, tOi3);

    data = zeros(nTrials*2, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    % data(nTrials*2+1:3*nTrials, :, :) = tempSig3;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];%, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN3(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(expv);
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

% Grand matrix concatenation (PEV) Omission Identity (A?B)

gkernel = ones(1, 1, 100); 

tOi1 = 1:4500;
tOi2 = 1:4500;
% tOi3 = 1:4500;
tN = length(tOi1);
icond1 = aaxb;
icond2 = bbxa;
% icond3 = rrrx;
nTrials = 100;

gmatrixN4 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, icond1};
    ncnt = size(tempSig1, 2);
    tcnt = size(tempSig1, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = tempSig1(iTrials, :, tOi1);

    tempSig2 = sspkData{ik, icond2};
    tcnt = size(tempSig2, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = tempSig2(iTrials, :, tOi2);

    % tempSig3 = sspkData{ik, icond3};
    % tcnt = size(tempSig3, 1);
    % iTrials = mod(randperm(nTrials), tcnt) + 1;
    % tempSig3 = tempSig3(iTrials, :, tOi3);

    data = zeros(nTrials*2, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    % data(nTrials*2+1:3*nTrials, :, :) = tempSig3;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];%, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN4(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(expv);
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Grand matrix concatenation (PEV) Stim Identity (A?B)

gkernel = ones(1, 1, 20); 

tOi1 = 501:1100;
tOi2 = 501:1100;
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

% tOi1 = 501:1500; %Azzz
tOi2 = 1531:2030; %zAzz
tOi3 = 2561:3060; %zzAz
tOi4 = 3591:4090; %zzzA

tN = length(tOi2);
nTrials = 100;

gmatrix2 = zeros(1, tN);
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, axab};
    tcnt = size(tempSig1, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig1 = tempSig1(iTrials, :, tOi2);
    ncnt = size(tempSig1, 2);

    tempSig2 = sspkData{ik, aaxb};
    tcnt = size(tempSig2, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig2 = tempSig2(iTrials, :, tOi3);

    tempSig3 = sspkData{ik, bbbx};
    tcnt = size(tempSig3, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig3 = tempSig3(iTrials, :, tOi4);

    tempSig4 = sspkData{ik, bxba};
    tcnt = size(tempSig4, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig4 = tempSig4(iTrials, :, tOi2);

    tempSig5 = sspkData{ik, bbxa};
    tcnt = size(tempSig5, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig5 = tempSig5(iTrials, :, tOi3);

    tempSig6 = sspkData{ik, aaax};
    tcnt = size(tempSig6, 1);
    iTrials = mod(randperm(nTrials), tcnt) + 1;
    tempSig6 = tempSig6(iTrials, :, tOi4);

    data = zeros(nTrials*6, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    data(2*nTrials+1:3*nTrials, :, :) = tempSig3;
    data(3*nTrials+1:4*nTrials, :, :) = tempSig4;
    data(4*nTrials+1:5*nTrials, :, :) = tempSig5;
    data(5*nTrials+1:6*nTrials, :, :) = tempSig6;

    groupIDs = [ones(1, 3*nTrials), ones(1, 3*nTrials)*2];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1, [1, 2], 1);
    gmatrix2(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(expv.*(p < 0.01));
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Grand matrix concatenation (iFR)

tOi1 = 1531:2030;
tOi2 = 1:500;
tN = length(tOi1);
icond1 = rrrr;
icond2 = rxrr;
nTrials = 100;

gmatrix3 = zeros(1, tN);
gmatrix4 = zeros(1, tN);
gmatrix5 = zeros(1, length(tOi2));
neuronCnt = 0;

for ik = 1:Nfiles

    tempSig1 = sspkData{ik, icond1};
    tempSig1 = mean(tempSig1(:, :, tOi1), 1);
   
    tempSig2 = sspkData{ik, icond2};
    ncnt = size(tempSig2, 2);
    tempSig2 = mean(tempSig2(:, :, tOi1), 1);
    
    tempSig3 = sspkData{ik, icond1};
    tempSig3 = mean(tempSig3(:, :, tOi2), 1);

    disp(ik);
    gmatrix3(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(tempSig1);
    gmatrix4(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(tempSig2);
    gmatrix5(neuronCnt+1:neuronCnt+ncnt, :) = squeeze(tempSig3);
    neuronCnt = neuronCnt + ncnt;
    disp(neuronCnt);

end

%% Neuron iFR plot (A)

icond1 = 1;
neuronID = 22;
fileID = 30;
neuronIDs = sum(fileIDs < fileID) + neuronID;

for nID = 1:length(neuronIDs)
    if nID == 1
        temp_sigx = squeeze(mean(sspkData{fileIDs(neuronIDs(nID)), icond1}(:, infileIDs(neuronIDs(nID)), :), 1));
        temp_sig1 = zeros([1, size(temp_sigx)]);
        temp_sig1(1, :, :) = temp_sigx;
    else
        temp_sig1(nID, :) = squeeze(mean(sspkData{fileIDs(neuronIDs(nID)), icond1}(:, infileIDs(neuronIDs(nID)), :), 1));
    end
end

figure;
imagesc(temp_sig1);

%% Single neuron iFR (N)

nID = 4099; % Grand neuron ID

icond1 = 1;
icond2 = 2;
icond3 = 3;
icond4 = 4;
kW = 500;
tN = 1000;%length(temp_sig1);

figure;

for icond = 1:4
    
    temp_sigx = sspkData{fileIDs(nID), icond}(:, infileIDs(nID), :);
    temp_sig1 = squeeze(mean(temp_sigx, 1));
    plot(tN*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", condNames(icond), "LineWidth", 2);
    hold("on");

end

xline(500, "HandleVisibility", "off");
xline(1530, "HandleVisibility", "off");
xline(2560, "HandleVisibility", "off");
xline(3590, "HandleVisibility", "off");

legend();
xlabel("Time(ms)");ylabel("FR(Spk/s)");
sgtitle("Neuron no." + num2str(nID) + " > " + areaList(areaIDs(nID)));

%% Single neuron PEV in time (N)

nID = 4099; % Grand neuron ID

kW = 10;

figure;

temp_sigx = gmatrixN2(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AAABvsBBBA)", "LineWidth", 2);
hold("on");

temp_sigx = gmatrixN1(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AAAXvsBBBX)", "LineWidth", 2);

temp_sigx = gmatrixN3(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AXABvsBXBA)", "LineWidth", 2);

temp_sigx = gmatrixN4(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AAXBvsBBXA)", "LineWidth", 2);


xline(500, "HandleVisibility", "off");
xline(1530, "HandleVisibility", "off");
xline(2560, "HandleVisibility", "off");
xline(3590, "HandleVisibility", "off");

legend();
xlabel("Time(ms)");ylabel("PEV%");
sgtitle("Neuron no." + num2str(nID) + " > " + areaList(areaIDs(nID)));
% >>> For the above neuron can you now plots its PEV time course?
% gmatrix0 >
% Same x-axis
% on y-axis, plot stimPEV (AAAB vs BBBA)
% and perhaps seperate PEV plots for stim identity:
% PEV for AXAB vs. BXBA vs. RXRR
% PEV for AAXB vs. BBXA vs. RRXR
% PEV for AAAX vs. BBBX vs. RRRX
% those 3 PEV plots could all be on the same plot
% maybe even together with the stim identity PEV
% use the same x axis (see how PEV changes over the trial) <<<

%%

mspkx = sspkData{1, 1};
mspkxo = mspkx;
mspkxo = smoothdata2(mspkxo, "gaussian", {1, 20});
mspkx = smoothdata2(mspkx, "gaussian", {1, 40});

kg = 12;
[idxs3, c1, sm1, d1] = kmeans(mspkx, kg, "Distance", "sqeuclidean");

figure;
sgtitle("All neurons>XA?XB>PEV>Kmeans");

for ik = 1:kg
    subplot(kg, 1, ik);
    plot(100*mean(mspkxo(idxs3 == ik, :)));
    xline(500, "Color", "#DD0000");
    xline(1530, "Color", "#DD0000");
    xline(2560, "Color", "#DD0000");
    xline(3590, "Color", "#DD0000");
    cntx = sum(idxs3 == ik);
    ylabel("PEV%");
    title("Kgroup" + num2str(ik) + ", nNeuron = " + num2str(cntx));
end

%%

% --- Script to Plot a Random Image in Multiple Tabs ---
% --- Create Multi-Tab Figure ---

fig = uifigure('Name', 'Multi-Tab Random Image Display', 'Position', [100 100 1700 1100]);
tabGroup = uitabgroup(fig, 'Position', [20 20 fig.Position(3)-40 fig.Position(4)-40]);

for i = 1:numberOfTabs
    % Create a new tab within the tab group
    tab = uitab(tabGroup, 'Title', ['Tab ' num2str(i)]);
    ax = uiaxes(tab);
    subplot(ax, 1, 1, 1);
    scatter(xe1(idxs), ye1(idxs), 1, color_t(iA, :), "filled", DisplayName="name");
    hold("on");

    title(ax, ['Random 100x100 Image - Plot ' num2str(i)]);

end

% --- End of Script ---

%%
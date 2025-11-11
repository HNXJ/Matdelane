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
sspkDataCleanTrials = cell(Nfiles, 12);

%% LFP unifier

% areax = "lfp";
% lfppath = "lfpSet\";
% lfpfiles = {dir(lfppath).name};
% lfpfiles = lfpfiles(contains(lfpfiles, areax));
% 
% Nfiles = length(lfpfiles);
% lfpData = cell(Nfiles, 1);

%% Load spk

% spkData = cell(Nfiles, 1);
% 
% for ik = 1:Nfiles
% 
%     tsetx = load(spkpath + spkfiles{ik});
%     tsetf = strsplit(spkfiles{ik}, "_");
%     spkArea{ik} = tsetf{2};
%     spkData{ik} = tsetx.xset;
% 
%     fprintf(num2str(ik) + "-");
% 
% end

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

% for ik = 1:Nfiles
% 
%     tsetx = load(lfppath + lfpfiles{ik});
%     lfpData{ik} = tsetx.xset;
%     fprintf(num2str(ik) + "-");
% 
% end

%% Find spike-less trials

for ik = 1:Nfiles

    ncnt = size(sspkData{ik, 1}, 2);
    tcnt = size(sspkData{ik, 1}, 3) / 1000;
    nafr = 1;

    for kk = 1:12

        trcnt = size(sspkData{ik, kk}, 1);
        sspkDataCleanTrials{ik, kk} = zeros(trcnt, ncnt);

    end

    for jk = 1:ncnt
        
        trnthsl = zeros(1, 12);
        trnthsu = zeros(1, 12);

        for kk = 1:12

            trcnt = size(sspkData{ik, kk}, 1);
            trn = squeeze(sum(sspkData{ik, kk}(:, jk, :), 3)) / tcnt;
            trnthsl(kk) = mean(trn(trn > median(trn))) - std(trn(trn > median(trn)));
            trnthsu(kk) = mean(trn(trn > median(trn))) + std(trn(trn > median(trn)));

        end

        trnthl = max(mean(trnthsl), 4.0);
        trnthu = max(mean(trnthsu), 6.0);

        for kk = 1:12

            trn = squeeze(sum(sspkData{ik, kk}(:, jk, :), 3)) / tcnt;
            trnthidx = (trn > trnthl) & (trn < trnthu);
            sspkDataCleanTrials{ik, kk}(trnthidx, jk) = 1;

        end

    end

    disp(ik);

end

%% Labels

areaList = ["V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"];
condNames = ["aaab", "axab", "aaxb", "aaax", ...
             "bbba", "bxba", "bbxa", "bbbx", ...
             "rrrr", "rxrr", "rrxr", "rrrx"];

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

neuronCnt = 0;

areaIDs = zeros(1, 1);
fileIDs = zeros(1, 1);
infileIDs = zeros(1, 1);
areaIDset = cell(1, 1);

for ik = 1:Nfiles
        
    disp(spkArea{ik});
    areaIDtemp = find(strcmpi(spkArea{ik}, areaList));
    ncnt = size(sspkData{ik, 1}, 2);
    areaIDset{ik} = neuronCnt+1:neuronCnt+ncnt;

    areaIDs(areaIDset{ik}) = areaIDtemp;
    fileIDs(areaIDset{ik}) = ik;
    infileIDs(areaIDset{ik}) = 1:ncnt;
    neuronCnt = neuronCnt + ncnt;

end

disp(neuronCnt);

%% Smoothing parameters

kW = 100;
gkernel = zeros(1, 1, kW);
gkernel(1, 1, :) = gausswin(kW);
pThresh = 0.01;

%% Grand matrix concatenation (IFR)

tOi1 = 1:4750;

tN = length(tOi1);

nTrials = 200;
gmatrixFR1 = zeros(12, neuronCnt, tN);

for ik = 1:Nfiles

    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials, ncnt, tN);

    for icond = 1:12

        for jk = 1:ncnt
    
            ltrials1 = find(sspkDataCleanTrials{ik, icond}(:, jk) == 1);
    
            tcnt1 = length(ltrials1);
    
            if tcnt1 > 10
    
                iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
                tempSig1(:, jk, :) = sspkData{ik, icond}(iTrials, jk, tOi1);
    
            end
    
        end
    
        gmatrixFR1(icond, areaIDset{ik}, :) = squeeze(mean(tempSig1, 1));

    end

    disp(ik);

end

%% Grand matrix concatenation (PEV) Stim Identity (A?B)

tOi1 = 1:4750;

tN = length(tOi1);
icond1 = aaab;
icond2 = bbba;
icond3 = rrrr;

nTrials = 200;
gmatrixN1 = zeros(neuronCnt, tN);

for ik = 1:Nfiles

    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials, ncnt, tN);
    tempSig2 = zeros(nTrials, ncnt, tN);
    tempSig3 = zeros(nTrials, ncnt, tN);

    for jk = 1:ncnt

        ltrials1 = find(sspkDataCleanTrials{ik, icond1}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, icond2}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, icond3}(:, jk) == 1);

        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);

        if tcnt1 > 10 & tcnt2 > 10 & tcnt3 > 10

            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSig1(:, jk, :) = sspkData{ik, icond1}(iTrials, jk, tOi1);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSig2(:, jk, :) = sspkData{ik, icond2}(iTrials, jk, tOi1);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSig3(:, jk, :) = sspkData{ik, icond3}(iTrials, jk, tOi1);

        end

    end

    data = zeros(nTrials*2, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    % data(2*nTrials+1:3*nTrials, :, :) = tempSig3;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];%, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN1(areaIDset{ik}, :) = squeeze(expv.*(p < pThresh));
    disp(ik);

end

% Grand matrix concatenation (PEV) 2nd Omission Identity (A?B)

tOi1 = 1:4750;

tN = length(tOi1);
icond1 = axab;
icond2 = axab;
icond3 = rxrr;

nTrials = 100;
gmatrixN2 = zeros(neuronCnt, tN);

for ik = 1:Nfiles

    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials, ncnt, tN);
    tempSig2 = zeros(nTrials, ncnt, tN);
    tempSig3 = zeros(nTrials, ncnt, tN);

    for jk = 1:ncnt

        ltrials1 = find(sspkDataCleanTrials{ik, icond1}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, icond2}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, icond3}(:, jk) == 1);

        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);

        if tcnt1 > 10 & tcnt2 > 10 & tcnt3 > 10

            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSig1(:, jk, :) = sspkData{ik, icond1}(iTrials, jk, tOi1);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSig2(:, jk, :) = sspkData{ik, icond2}(iTrials, jk, tOi1);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSig3(:, jk, :) = sspkData{ik, icond3}(iTrials, jk, tOi1);

        end

    end

    data = zeros(nTrials*3, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    data(2*nTrials+1:3*nTrials, :, :) = tempSig3;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN2(areaIDset{ik}, :) = squeeze(expv.*(p < pThresh));
    disp(ik);

end

% Grand matrix concatenation (PEV) 3rd Omission Identity (A?B)

tOi1 = 1:4750;

tN = length(tOi1);
icond1 = aaxb;
icond2 = aaxb;
icond3 = rrxr;

nTrials = 100;
gmatrixN3 = zeros(neuronCnt, tN);

for ik = 1:Nfiles

    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials, ncnt, tN);
    tempSig2 = zeros(nTrials, ncnt, tN);
    tempSig3 = zeros(nTrials, ncnt, tN);

    for jk = 1:ncnt

        ltrials1 = find(sspkDataCleanTrials{ik, icond1}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, icond2}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, icond3}(:, jk) == 1);

        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);

        if tcnt1 > 10 & tcnt2 > 10 & tcnt3 > 10

            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSig1(:, jk, :) = sspkData{ik, icond1}(iTrials, jk, tOi1);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSig2(:, jk, :) = sspkData{ik, icond2}(iTrials, jk, tOi1);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSig3(:, jk, :) = sspkData{ik, icond3}(iTrials, jk, tOi1);

        end

    end

    data = zeros(nTrials*3, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    data(2*nTrials+1:3*nTrials, :, :) = tempSig3;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN3(areaIDset{ik}, :) = squeeze(expv.*(p < pThresh));
    disp(ik);

end

% Grand matrix concatenation (PEV) 4th Omission Identity (A?B)

tOi1 = 1:4750;

tN = length(tOi1);
icond1 = aaax;
icond2 = aaax;
icond3 = rrrx;

nTrials = 100;
gmatrixN4 = zeros(neuronCnt, tN);

for ik = 1:Nfiles

    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials, ncnt, tN);
    tempSig2 = zeros(nTrials, ncnt, tN);
    tempSig3 = zeros(nTrials, ncnt, tN);

    for jk = 1:ncnt

        ltrials1 = find(sspkDataCleanTrials{ik, icond1}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, icond2}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, icond3}(:, jk) == 1);

        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);

        if tcnt1 > 10 & tcnt2 > 10 & tcnt3 > 10

            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSig1(:, jk, :) = sspkData{ik, icond1}(iTrials, jk, tOi1);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSig2(:, jk, :) = sspkData{ik, icond2}(iTrials, jk, tOi1);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSig3(:, jk, :) = sspkData{ik, icond3}(iTrials, jk, tOi1);

        end

    end

    data = zeros(nTrials*3, ncnt, tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    data(2*nTrials+1:3*nTrials, :, :) = tempSig3;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1);
    gmatrixN4(areaIDset{ik}, :) = squeeze(expv.*(p < pThresh));
    disp(ik);

end

%% Grand matrix concatenation (PEV) Stim Identity (A?B)

tOi1 = 1:3500;

tN = length(tOi1);
icond1 = aaab;
icond2 = bbba;
nTrials = 100;

gmatrix1 = zeros(neuronCnt, tN);

for ik = 1:Nfiles

    nTrials = size(sspkData{ik, icond1}, 1);
    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials, ncnt, tN);
    tempSig2 = zeros(nTrials, ncnt, tN);

    for jk = 1:size(tempSig1, 2)

        ltrials1 = find(sspkDataCleanTrials{ik, icond1}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, icond2}(:, jk) == 1);
        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);

        if tcnt1 > 10 & tcnt2 > 10

            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSig1(:, jk, :) = sspkData{ik, icond1}(iTrials, jk, tOi1);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSig2(:, jk, :) = sspkData{ik, icond2}(iTrials, jk, tOi1);

        end

    end

    data = zeros(nTrials*2, size(tempSig2, 2), tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1, [1, 2]);
    gmatrix1(areaIDset{ik}, :) = squeeze(expv.*(p < pThresh));
    disp(ik);

end

%% Grand matrix concatenation (PEV) Omission Identity (X|A?X|B)

% tOi1 = 501:1050; %Azzz
tOi2 = 1501:2300; %ZAzz
tOi3 = 2531:3330; %zZAz
tOi4 = 3561:4360; %zzZA

tN = length(tOi2);

gmatrix2 = zeros(neuronCnt, tN);

for ik = 1:Nfiles

    nTrials = size(sspkData{ik, axab}, 1);
    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials*3, ncnt, tN);
    tempSig2 = zeros(nTrials*3, ncnt, tN);
        
    tOi2shuffle = tOi2(randperm(length(tOi2)));
    tOi3shuffle = tOi3(randperm(length(tOi3)));
    tOi4shuffle = tOi4(randperm(length(tOi4)));

    for jk = 1:size(tempSig1, 2)

        ltrials1 = find(sspkDataCleanTrials{ik, axab}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, aaxb}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, bbbx}(:, jk) == 1);
        ltrials4 = find(sspkDataCleanTrials{ik, bxba}(:, jk) == 1);
        ltrials5 = find(sspkDataCleanTrials{ik, bbxa}(:, jk) == 1);
        ltrials6 = find(sspkDataCleanTrials{ik, aaax}(:, jk) == 1);

        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);
        tcnt4 = length(ltrials4);
        tcnt5 = length(ltrials5);
        tcnt6 = length(ltrials6);

        if tcnt1 > 10 & tcnt2 > 10 & tcnt3 > 10 & tcnt4 > 10 & tcnt5 > 10 & tcnt6 > 10

            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSig1(1:nTrials, jk, :) = sspkData{ik, axab}(iTrials, jk, tOi2);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSig1(nTrials+1:2*nTrials, jk, :) = sspkData{ik, aaxb}(iTrials, jk, tOi3);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSig1(2*nTrials+1:3*nTrials, jk, :) = sspkData{ik, bbbx}(iTrials, jk, tOi4);

            iTrials = ltrials4(mod(1:nTrials, tcnt4) + 1);
            tempSig2(1:nTrials, jk, :) = sspkData{ik, bxba}(iTrials, jk, tOi2shuffle);
            iTrials = ltrials5(mod(1:nTrials, tcnt5) + 1);
            tempSig2(nTrials+1:2*nTrials, jk, :) = sspkData{ik, bbxa}(iTrials, jk, tOi3shuffle);
            iTrials = ltrials6(mod(1:nTrials, tcnt6) + 1);
            tempSig2(2*nTrials+1:3*nTrials, jk, :) = sspkData{ik, aaax}(iTrials, jk, tOi4shuffle);
        end

    end

    data = zeros(nTrials*6, size(tempSig1, 2), tN);
    data(1:3*nTrials, :, :) = tempSig1;
    data(3*nTrials+1:6*nTrials, :, :) = tempSig2;

    groupIDs = [ones(1, 3*nTrials), ones(1, 3*nTrials)*2];

    data = convn(data, gkernel, 'same');
    % data = squeeze(mean(data, 3));
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1, [1, 2]);
    gmatrix2(areaIDset{ik}, :) = squeeze(expv.*(p < pThresh));
    disp(ik);

end

%% Grand matrix concatenation (PEV) Omission Position (X|s2,s3,s4)

% tOi1 = 501:1500; %Azzz
tOi2 = 1531:2030; %zAzz
tOi3 = 2561:3060; %zzAz
tOi4 = 3591:4090; %zzzA

tN = length(tOi2);
nTrials = 100;

gmatrix6 = zeros(neuronCnt, tN);

for ik = 1:Nfiles

    ncnt = length(areaIDset{ik});
    tempSig1 = zeros(nTrials, ncnt, tN);

    for jk = 1:size(tempSig1, 2)
        ltrials = find(sspkDataCleanTrials{ik, axab}(:, jk) == 1);
        tcnt = length(ltrials);
        if tcnt > 10
            iTrials = ltrials(mod(randperm(nTrials), tcnt) + 1);
            tempSig1(:, jk, :) = sspkData{ik, axab}(iTrials, jk, tOi2);
        end
    end

    tempSig2 = zeros(nTrials, ncnt, tN);

    for jk = 1:size(tempSig2, 2)
        ltrials = find(sspkDataCleanTrials{ik, aaxb}(:, jk) == 1);
        tcnt = length(ltrials);
        if tcnt > 10
            iTrials = ltrials(mod(randperm(nTrials), tcnt) + 1);
            tempSig2(:, jk, :) = sspkData{ik, aaxb}(iTrials, jk, tOi3);
        end
    end

    tempSig3 = zeros(nTrials, ncnt, tN);

    for jk = 1:size(tempSig3, 2)
        ltrials = find(sspkDataCleanTrials{ik, bbbx}(:, jk) == 1);
        tcnt = length(ltrials);
        if tcnt > 10
            iTrials = ltrials(mod(randperm(nTrials), tcnt) + 1);
            tempSig3(:, jk, :) = sspkData{ik, bbbx}(iTrials, jk, tOi4);
        end
    end

    tempSig4 = zeros(nTrials, ncnt, tN);

    for jk = 1:size(tempSig4, 2)
        ltrials = find(sspkDataCleanTrials{ik, bxba}(:, jk) == 1);
        tcnt = length(ltrials);
        if tcnt > 10
            iTrials = ltrials(mod(randperm(nTrials), tcnt) + 1);
            tempSig4(:, jk, :) = sspkData{ik, bxba}(iTrials, jk, tOi2);
        end
    end

    tempSig5 = zeros(nTrials, ncnt, tN);

    for jk = 1:size(tempSig5, 2)
        ltrials = find(sspkDataCleanTrials{ik, bbxa}(:, jk) == 1);
        tcnt = length(ltrials);
        if tcnt > 10
            iTrials = ltrials(mod(randperm(nTrials), tcnt) + 1);
            tempSig5(:, jk, :) = sspkData{ik, bbxa}(iTrials, jk, tOi3);
        end
    end

    tempSig6 = zeros(nTrials, ncnt, tN);

    for jk = 1:size(tempSig6, 2)
        ltrials = find(sspkDataCleanTrials{ik, aaax}(:, jk) == 1);
        tcnt = length(ltrials);
        if tcnt > 10
            iTrials = ltrials(mod(randperm(nTrials), tcnt) + 1);
            tempSig6(:, jk, :) = sspkData{ik, aaax}(iTrials, jk, tOi4);
        end
    end

    data = zeros(nTrials*6, size(tempSig1, 2), tN);
    data(1:nTrials, :, :) = tempSig1;
    data(nTrials+1:2*nTrials, :, :) = tempSig2;
    data(2*nTrials+1:3*nTrials, :, :) = tempSig3;

    data(3*nTrials+1:4*nTrials, :, :) = tempSig4;
    data(4*nTrials+1:5*nTrials, :, :) = tempSig5;
    data(5*nTrials+1:6*nTrials, :, :) = tempSig6;
    groupIDs = [ones(1, nTrials), ones(1, nTrials)*2, ones(1, nTrials)*3 ...
        , ones(1, nTrials)*1, ones(1, nTrials)*2, ones(1, nTrials)*3];

    data = convn(data, gkernel, 'same');
    [expv, n, mu, p, F] = jPEV(data, groupIDs, 1, [1, 2, 3], 1);
    gmatrix6(areaIDset{ik}, :) = squeeze(expv.*(p < pThresh));
    disp(ik);

end

%% Grand matrix concatenation (iFR)

tOi1 = 501:1300; %Azzz
tOi2 = 1531:2330; %zAzz
tOi3 = 2561:3360; %zzAz
tOi4 = 3591:4390; %zzzA

tOib = 1:400;
tOis = [501:1000, 1531:2030, 2561:3060, 3591:4090];

tN = length(tOis);
nTrials = 100;
icond1 = rrrr;

gmatrix3 = zeros(neuronCnt, tN);
gmatrix4 = zeros(neuronCnt, length(tOi1));
gmatrix5 = zeros(neuronCnt, length(tOib));

for ik = 1:Nfiles

    ncnt = length(areaIDset{ik});
    tempSigG = zeros(nTrials, ncnt, size(sspkData{ik, icond1}, 3));

    for jk = 1:size(tempSigG, 2)
        ltrials = find(sspkDataCleanTrials{ik, icond1}(:, jk) == 1);
        tcnt = length(ltrials);
        if tcnt > 10
            iTrials = ltrials(mod(1:nTrials, tcnt) + 1);
            tempSigG(:, jk, :) = sspkData{ik, icond1}(iTrials, jk, :);
        end
    end

    tempSigx1 = zeros(nTrials*3, ncnt, size(sspkData{ik, icond1}, 3));

    for jk = 1:size(tempSigx1, 2)
        ltrials1 = find(sspkDataCleanTrials{ik, axab}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, bxba}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, rxrr}(:, jk) == 1);
        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);
        if  tcnt1 > 5 & tcnt2 > 5 & tcnt3 > 5
            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSigx1(1:nTrials, jk, :) = sspkData{ik, axab}(iTrials, jk, :);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSigx1(nTrials*1+1:2*nTrials, jk, :) = sspkData{ik, bxba}(iTrials, jk, :);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSigx1(nTrials*2+1:3*nTrials, jk, :) = sspkData{ik, rxrr}(iTrials, jk, :);
        end
    end
        
    tempSigx2 = zeros(nTrials*3, ncnt, size(sspkData{ik, icond1}, 3));

    for jk = 1:size(tempSigx1, 2)
        ltrials1 = find(sspkDataCleanTrials{ik, aaxb}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, bbxa}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, rrxr}(:, jk) == 1);
        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);
        if tcnt1 > 5 & tcnt2 > 5 & tcnt3 > 5
            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSigx2(1:nTrials, jk, :) = sspkData{ik, aaxb}(iTrials, jk, :);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSigx2(nTrials*1+1:2*nTrials, jk, :) = sspkData{ik, bbxa}(iTrials, jk, :);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSigx2(nTrials*2+1:3*nTrials, jk, :) = sspkData{ik, rrxr}(iTrials, jk, :);
        end
    end

    tempSigx3 = zeros(nTrials*3, ncnt, size(sspkData{ik, icond1}, 3));

    for jk = 1:size(tempSigx1, 2)
        ltrials1 = find(sspkDataCleanTrials{ik, aaax}(:, jk) == 1);
        ltrials2 = find(sspkDataCleanTrials{ik, bbbx}(:, jk) == 1);
        ltrials3 = find(sspkDataCleanTrials{ik, rrrx}(:, jk) == 1);
        tcnt1 = length(ltrials1);
        tcnt2 = length(ltrials2);
        tcnt3 = length(ltrials3);
        if tcnt1 > 5 & tcnt2 > 5 & tcnt3 > 5
            iTrials = ltrials1(mod(1:nTrials, tcnt1) + 1);
            tempSigx3(1:nTrials, jk, :) = sspkData{ik, aaax}(iTrials, jk, :);
            iTrials = ltrials2(mod(1:nTrials, tcnt2) + 1);
            tempSigx3(nTrials*1+1:2*nTrials, jk, :) = sspkData{ik, bbbx}(iTrials, jk, :);
            iTrials = ltrials3(mod(1:nTrials, tcnt3) + 1);
            tempSigx3(nTrials*2+1:3*nTrials, jk, :) = sspkData{ik, rrrx}(iTrials, jk, :);
        end
    end

    tempSig1 = mean(tempSigG(:, :, tOis), 1);
    tempSig2 = (mean(tempSigx1(:, :, tOi2), 1) + mean(tempSigx2(:, :, tOi3), 1) + mean(tempSigx3(:, :, tOi4), 1))/3;
    tempSig3 = mean(tempSigG(:, :, tOib), 1);

    gmatrix3(areaIDset{ik}, :) = squeeze(tempSig1);
    gmatrix4(areaIDset{ik}, :) = squeeze(tempSig2);
    gmatrix5(areaIDset{ik}, :) = squeeze(tempSig3);
    disp(ik);

end

%% Neuron iFR plot (A)

icond1 = 1;
neuronID = 1:100;
fileID = 31;
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
xlabel("time(ms)");
ylabel("Neuron no.");
title("Rastrogram");

%% Single neuron iFR (N)

nID = 236; % Grand neuron ID 4099(FST) | 3461(FEF)

kW = 250;
kX = 1;
tN = 1000; % length(temp_sig1);
timevec = linspace(-500, 4250, 4750);

figure;
condOi = [4, 8];
ncondOi = length(condOi);

for ik = 1:ncondOi
    
    icond = condOi(ik);
    temp_sigx = squeeze(sspkData{fileIDs(nID), icond}(:, infileIDs(nID), :));
    temp_sigx = tN*smoothdata2(temp_sigx, "gaussian", {1, kW});

    temp_sig1 = squeeze(mean(temp_sigx, 1));
    temp_sig2 = squeeze(std(temp_sigx, 1) / sqrt(size(temp_sigx, 1)));
    % temp_sig1 = tN*smoothdata(temp_sig1, "gaussian", kW);
    % temp_sig2 = tN*smoothdata(temp_sig2, "gaussian", 1);

    xpe1 = timevec;
    ype1 = temp_sig1 + kX*temp_sig2;
    xpe2 = timevec(end:-1:1);
    ype2 = temp_sig1(end:-1:1) - kX*temp_sig2(end:-1:1);

    plot(timevec, temp_sig1, "DisplayName", condNames(icond) + "+-" + num2str(kX) + "SEM", "LineWidth", 2, "Color", color_t(icond, :));
    patch([xpe1, xpe2], [ype1, ype2], color_t(icond, :), "FaceAlpha", 1.0, "HandleVisibility", "off", "EdgeColor", "none");
    hold("on");

end

xline(0, "HandleVisibility", "off");
xline(1030, "HandleVisibility", "off");
xline(2060, "HandleVisibility", "off");
xline(3090, "HandleVisibility", "off");

xlim([-750, 4500]);
legend();
xlabel("Time(ms)");ylabel("FR(Spk/s)");
sgtitle("Neuron no." + num2str(nID) + " > " + areaList(areaIDs(nID)));

%% Single neuron PEV in time (N)

nID = 3007; % Grand neuron ID

kW = 100;

figure;

temp_sigx = gmatrixN1(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AAABvsBBBA)", "LineWidth", 2);
hold("on");

temp_sigx = gmatrixN2(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AXAXvsBXBA)", "LineWidth", 2);

temp_sigx = gmatrixN3(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AAXBvsBBXA)", "LineWidth", 2);

temp_sigx = gmatrixN4(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
plot(100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", "PEV(AAAXvsBBBX)", "LineWidth", 2);

xline(500, "HandleVisibility", "off");
xline(1530, "HandleVisibility", "off");
xline(2560, "HandleVisibility", "off");
xline(3590, "HandleVisibility", "off");

legend();
xlabel("Time(ms)");ylabel("PEV%");
sgtitle("Neuron no." + num2str(nID) + " > " + areaList(areaIDs(nID)) + " > smoothW = " + num2str(kW));

%% Single neuron rastrogram (N)

nID = 34; % Grand neuron ID % 1004/1269

icond1 = 1;
icond2 = 2;
icond3 = 3;
icond4 = 4;

kW = 1;
kX = 1;
tN = 1*kW; % length(temp_sig1);
timevec = linspace(-500, 4250, 4750);

figure;
condOi =[1, 4, 8, 12];
ncondOi = length(condOi);

for ik = 1:ncondOi
    
    icond = condOi(ik);
    temp_sigxc = find(squeeze(sspkDataCleanTrials{fileIDs(nID), icond}(:, infileIDs(nID))) == 1);
    temp_sigx = squeeze(sspkData{fileIDs(nID), icond}(temp_sigxc, infileIDs(nID), :));
    temp_sigx = tN*smoothdata2(temp_sigx, "gaussian", {1, kW});

    subplot(ceil(ncondOi/2), 2, ik);
    imagesc(temp_sigx, "XData", timevec);hold("on");
    title(condNames(icond));
    ylabel("Trial no.");

    xline(0, "HandleVisibility", "off", "Color", [1 1 1]);
    xline(1030, "HandleVisibility", "off", "Color", [0 1 1]);
    xline(2060, "HandleVisibility", "off", "Color", [0 1 1]);
    xline(3090, "HandleVisibility", "off", "Color", [0 1 1]);

    % cb = colorbar();
    % ylabel(cb, "Spk/s");

end

xlabel("Time(ms)");
sgtitle("Neuron no." + num2str(nID) + " > " + areaList(areaIDs(nID)));

%% Single neuron PEV in time (N) with iFR

nID = 34; % Grand neuron ID 4099(FST) | 3461(FEF) | 1106//4094 | 3602 | MST 3007

kW = 400;
kX = 2;
tN = 1000; % length(temp_sig1);
timevec = linspace(-500, 4250, 4750);

figure;
condOi = [2, 3, 6, 7, 10, 11];
ncondOi = length(condOi);

for ik = 1:ncondOi
    
    icond = condOi(ik);

    temp_sigxc = find(squeeze(sspkDataCleanTrials{fileIDs(nID), icond}(:, infileIDs(nID))) == 1);
    temp_sigx = squeeze(sspkData{fileIDs(nID), icond}(temp_sigxc, infileIDs(nID), :));
    temp_sigx = tN*smoothdata2(temp_sigx, "gaussian", {1, kW});

    temp_sig1 = squeeze(mean(temp_sigx, 1));
    temp_sig2 = squeeze(std(temp_sigx, 1) / sqrt(size(temp_sigx, 1)));
    % temp_sig1 = tN*smoothdata(temp_sig1, "gaussian", kW);
    % temp_sig2 = tN*smoothdata(temp_sig2, "gaussian", 1);

    xpe1 = timevec;
    ype1 = temp_sig1 + kX*temp_sig2;
    xpe2 = timevec(end:-1:1);
    ype2 = temp_sig1(end:-1:1) - kX*temp_sig2(end:-1:1);

    % plot(timevec, temp_sig1, "LineWidth", 2, "Color", color_t(icond, :));
    patch([xpe1, xpe2], [ype1, ype2], color_t(icond, :), "FaceAlpha", 0.2, "EdgeColor", "none", "DisplayName", condNames(icond) + "+-" + num2str(kX) + "SEM");
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
dispnametemp1 = "PEV(AXABvsBXBA)";
temp_sigx = gmatrixN2(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
yyaxis("right");

% kW = 1;
plot(timevec, 100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", dispnametemp1, "LineWidth", 2, "color", [.5 0 .5]);
yline(0, "LineStyle", "--");

dispnametemp1 = "PEV(AAXBvsBBXA)";
temp_sigx = gmatrixN3(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
yyaxis("right");

% kW = 1;
plot(timevec, 100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", dispnametemp1, "LineWidth", 2, "color", [.5 0 .5]);
yline(0, "LineStyle", "--");

dispnametemp1 = "PEV(AAAXvsBBBX)";
temp_sigx = gmatrixN4(nID, :);
temp_sig1 = squeeze(mean(temp_sigx, 1));
yyaxis("right");

% kW = 1;
plot(timevec, 100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", dispnametemp1, "LineWidth", 2, "color", [.5 0 .5]);
yline(0, "LineStyle", "--");

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ylabel("PEV(%)", "Color", [1 0 0]);
sgtitle("Neuron no." + num2str(nID) + " > " + areaList(areaIDs(nID)));

%% Clustering (Kmeans)

mspkx = sspkData{35, 10};
mspkx = squeeze(mean(mspkx, 1));
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

%% Coloring

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

%% Group PEVs in time (N) with iFRs

% nIDgroup = pfcidxng1;
nIDgroup = fefidxng1;

kW = 200;
kX = 1;
tN = 1000; % length(temp_sig1);
timevec = linspace(-500, 4250, 4750);

figure;
condOi = [2, 6, 10];
ncondOi = length(condOi);

for ik = 1:ncondOi
    
    icond = condOi(ik);
    temp_sig1 = zeros(size(timevec));
    temp_sig2 = zeros(size(timevec));

    for nIDx = nIDgroup

        temp_sigxc = find(squeeze(sspkDataCleanTrials{fileIDs(nIDx), icond}(:, infileIDs(nIDx))) == 1);
        temp_sigx = squeeze(sspkData{fileIDs(nIDx), icond}(temp_sigxc, infileIDs(nIDx), :));
        temp_sigx = tN*smoothdata2(temp_sigx, "gaussian", {1, kW});
    
        temp_sig1 = temp_sig1 + squeeze(mean(temp_sigx, 1))/length(nIDgroup);
        temp_sig2 = temp_sig2 + squeeze(std(temp_sigx, 1) / sqrt(size(temp_sigx, 1)))/length(nIDgroup);

    end

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
temp_sig1 = squeeze(mean(temp_sigx, 1));
temp_sig1 = temp_sig1 - min(temp_sig1);
yyaxis("right");

% kW = 1;
plot(timevec, 100*smoothdata(temp_sig1, "gaussian", kW), "DisplayName", dispnametemp1, "LineWidth", 2, "color", [.5 0 .5]);
yline(0, "LineStyle", "--");

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ylabel("PEV(%)", "Color", [1 0 0]);
sgtitle("N = " + num2str(length(nIDgroup)) + " > " + areaList(areaIDs(nIDgroup(1))));

%%
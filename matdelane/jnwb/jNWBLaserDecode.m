function [starttime_ms, conditions] =  jNWBLaserDecode(nwb)
%laser signal and timestamps
laser1 = nwb.acquisition.get('laser_1_tracking').timeseries.get('laser_1_tracking_data');
lasersignal = laser1.data.load();
lasersmooth = smoothdata(lasersignal, 'gaussian', 3);
lasertimestamps = laser1.timestamps.load();

%detect when laser turned on (oscillates around ~0.02 between trials)
threshold = 0.1;
starttime = ((find(diff([0; lasersignal > threshold]) == 1))-1) / 1000;
starttime_ms = floor(starttime*1000);
stoptimes = (find(diff([0; lasersignal > threshold]) == -1)) / 1000;
stoptimes_ms = floor(stoptimes*1000);
trialnum = 1:length(starttime);

% check trial durations, which should be 2s
trlduration = stoptimes_ms-starttime_ms;
shorttrl = find(trlduration < 2000);

% cut out short trials, this at least works for 240711,240712
starttime_ms(shorttrl) = [];
stoptimes_ms(shorttrl) = [];
starttime(shorttrl) = [];
stoptimes(shorttrl) = [];
trialnum(shorttrl) = [];
numtrials = length(starttime);
% but now you have stop times that are all offset because
% the signal stays at ~0.7 for ~15-30ms before returning to baseline(off)
stoptimes = starttime+2;
stoptimes_ms = floor(stoptimes*1000);


% amplitudes
decodedamptrl = zeros(1, length(starttime));
for i = 1:numtrials
    pulse = lasersignal(starttime_ms(i):stoptimes_ms(i));
    decodedamptrl(i) = round(max(pulse),1);
end
decodedamplist = unique(decodedamptrl);

% fft to decode frequency for each trial
P2_all = cell(numtrials, 1);
P1_all = cell(numtrials, 1);
f_all = cell(numtrials, 1);
for i = 1:numtrials
    pulse = lasersmooth(starttime_ms(i):stoptimes_ms(i));
    pulse = pulse - mean(pulse);
    trlinterval = length(pulse);
    trlfft = fft(pulse);
    P2 = abs(trlfft/trlinterval);
    P1 = P2(1:floor(trlinterval/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (0:(trlinterval/2))*(1000/trlinterval);
    P2_all{i} = P2;
    P1_all{i} = P1;
    f_all{i} = f;
end

decodedfreqtrl = zeros(numtrials, 1);
for i = 1:numtrials
    [peakValue, peakIndex] = max(P1_all{i});
    decodedfreqtrl(i) = round(f_all{i}(peakIndex) * 2) / 2; % Round to nearest 0.5
end
decodedfreqlist = unique(decodedfreqtrl);

% write conditions array
% check how many conditions this session has
nwbconds = nwb.intervals.get('passive_opto').vectordata.get("task_condition_number").data(:);
nwbcondlist = unique(nwbconds);
% and change this if needed
numconds = length(nwbcondlist);
switch numconds
    case 24
        freq_conds = [5, 11, 21, 35, 37.5, 40, 75, 115];
        amp_conds = [2.3, 3.8, 5.3];
        ctrl_amp = NaN;
    case 25
        freq_conds = [5, 11, 21, 35, 37.5, 40, 75, 115]; 
        amp_conds = [3.7, 4.9, 6.1];
        ctrl_amp = 0.7;
    case 21
        freq_conds = [2, 5, 11, 21, 31, 38, 41, 75, 115, 140];
        amp_conds = [3.7, 6.1];
        ctrl_amp = 0.7;
end

conditions = zeros(size(trialnum));
for i = 1:length(trialnum)
    if  ~isnan(ctrl_amp) && abs(decodedamptrl(i) - ctrl_amp) < 0.2
        conditions(i) = length(freq_conds) * length(amp_conds) + 1;
    else
        [~, freq_idx] = min(abs(freq_conds - decodedfreqtrl(i)));
        [~, amp_idx] = min(abs(amp_conds - decodedamptrl(i)));
        conditions(i) = (amp_idx - 1) * length(freq_conds) + freq_idx;
    end
end

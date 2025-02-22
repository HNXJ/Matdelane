function jNWBLaserTrialDecode(nwb)
% nwbFile = 'sub-C31s_ses-240730.nwb';
% nwb = nwbRead(nwbFile);

%laser signal and timestamps
laser1 = nwb.acquisition.get('laser_1_tracking').timeseries.get('laser_1_tracking_data');
lasersignal = laser1.data.load();
lasersmooth = smoothdata(lasersignal, 'gaussian', 3);
lasertimestamps = laser1.timestamps.load();

%detect when laser turned on/off 
threshold = 0.1;
starttime = ((find(diff([0; lasersignal > threshold]) == 1))-1) / 1000;
starttime_ms = floor(starttime*1000);
stoptimes = (find(diff([0; lasersignal > threshold]) == -1)) / 1000;
stoptimes_ms = floor(stoptimes*1000);
trialnum = 1:length(starttime);
numtrials = length(starttime);

% check for trials where monkey broke fixation early
trldur = stoptimes_ms-starttime_ms;
shorttrl = find(trldur < 1999);
correcttrl = trldur >= 1999;
% exclude those trials
starttime_ms(shorttrl) = [];
stoptimes_ms(shorttrl) = [];
starttime(shorttrl) = [];
stoptimes(shorttrl) = [];
trialnum(shorttrl) = [];
numtrials = length(starttime);
trldur = stoptimes_ms-starttime_ms;

%round down trls where signal stays at 0.7 baseline for ~30ms at the end
stoptimes = starttime+2;
stoptimes_ms = floor(stoptimes*1000);

% amplitudes
decodedamptrlraw = zeros(1, length(starttime));
for i = 1:numtrials
    pulse = lasersignal(starttime_ms(i):stoptimes_ms(i));
    decodedamptrlraw(i) = round(max(pulse),1);
end
decodedamplistraw = unique(decodedamptrlraw);


% baseline corrected rounded amplitudes
decodedamptrl = zeros(1, length(starttime));
for i = 1:numtrials
    pulse = lasersignal(starttime_ms(i):stoptimes_ms(i));
    decodedamptrl(i) = round(max(pulse),1);
end
% subtract baseline
if any(decodedamptrl < 1)
    baseline_amp = min(decodedamptrl);
else
    baseline_amp = 0.7;
end
decodedamptrl = decodedamptrl - baseline_amp;
decodedamptrl = round(decodedamptrl);
% decodedamptrl = round(decodedamptrl / 0.1) * 0.1;
decodedamplist = unique(decodedamptrl);

% fft
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
% frequencies for all non-control trials. if monkey broke fixation before a
% full cycle of the laser, frequency for that trial is inaccurate
decodedfreqtrl = NaN(numtrials, 1);
for i = 1:numtrials
    if decodedamptrl(i) > 1 
        [peakValue, peakIndex] = max(P1_all{i});
        decodedfreqtrl(i) = round(f_all{i}(peakIndex));
    end
end
decodedfreqlist = unique(decodedfreqtrl(~isnan(decodedfreqtrl)));



% define condition sets
cond_sets = {
    struct('freqs', [5, 11, 20, 35, 37, 40, 75, 115], 'amps', [0, 3, 4, 5]), ...
    struct('freqs', [5, 11, 20, 35, 37, 40, 75, 115], 'amps', [2, 3, 5]), ...
    struct('freqs', [5, 11, 21, 32, 41, 71, 101, 140], 'amps', [0, 3, 4, 5]), ...
    struct('freqs', [1, 5, 11, 21, 31, 38, 40, 75, 105, 140], 'amps', [0, 3, 5]), ...
    struct('freqs', [2, 5, 11, 21, 31, 38, 40, 75, 105, 140], 'amps', [0, 3, 5])
};

% compare decoded freq and ampl to condition sets
selected_set = [];
for j = 1:length(cond_sets)
    if all(ismember(decodedfreqlist, cond_sets{j}.freqs)) && ...
       all(ismember(decodedamplist, cond_sets{j}.amps))
        selected_set = cond_sets{j};
        break;
    end
end

if isempty(selected_set)
    error('No matching condition set found for the decoded data.');
end

freq_conds = selected_set.freqs;
amp_conds = selected_set.amps;

% map conditions to trials
laservoltage = amp_conds(amp_conds > 0);
conditions = zeros(size(trialnum));
for i = 1:length(trialnum)
    if decodedamptrl(i) == 0
        conditions(i) = length(freq_conds) * length(laservoltage) + 1;
    else
        [~, freq_idx] = min(abs(freq_conds - decodedfreqtrl(i)));
        [~, amp_idx] = min(abs(laservoltage - decodedamptrl(i)));
        conditions(i) = (amp_idx - 1) * length(freq_conds) + freq_idx;
    end
end
conditionslist = unique(conditions);


% plot trials x conditions and verify 
figure(1);
tiledlayout('flow');
sgtitle(nwbFile);
t = 4000;
epoched_laser = zeros(numtrials, t);
for cond_idx = 11:21
    cond = conditionslist(cond_idx);
    curr_cond_idx = find(conditions == cond);
    numtrial_cond = length(curr_cond_idx);
    epoched_laser_cond = zeros(numtrial_cond, t);
    
    if cond == length(freq_conds) * length(laservoltage) + 1
        freq_label = 'Control';
        amp_label = '';
    else
        freq_idx = mod(cond - 1, length(freq_conds)) + 1;
        amp_idx = floor((cond - 1) / length(freq_conds)) + 1;
        freq_label = sprintf('%d Hz', freq_conds(freq_idx));
        amp_label = sprintf('%.1f V', laservoltage(amp_idx));
    end
    
    for trl = 1:numtrial_cond
        epoched_laser_cond(trl, :) = lasersignal(starttime_ms(curr_cond_idx(trl)) - 1000 : ...
                                                 starttime_ms(curr_cond_idx(trl)) + 2999);
    end

    figure(1);
    nexttile;
    imagesc(epoched_laser_cond);
    if strcmp(freq_label, 'Control')
        title(sprintf('%s, %d trials', freq_label, numtrial_cond));
    else
        title(sprintf('%s, %s, %d trials', freq_label, amp_label, numtrial_cond));
    end
end

filename = nwbFile;
filename = replace(filename, ".nwb", ".mat");
save(filename, "conditions", "conditionslist", "trialnum", "starttime","stoptimes",...
    "starttime_ms","stoptimes_ms","decodedamplist","decodedfreqlist",...
    "decodedamptrl","decodedfreqtrl", "correcttrl");

end
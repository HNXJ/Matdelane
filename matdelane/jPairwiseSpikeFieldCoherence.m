function SPKF = jPairwiseSpikeFieldCoherence(LFP, SPK, Fs, fRange)
    % Inputs:
    %   LFP: Local Field Potential data [Trial x Channel x Time]
    %   SPK: Spiking data [Trial x NeuronID x Time] (binary or spike times)
    %   Fs: Sampling frequency (in Hz)
    %   fRange: Frequency range for coherence [fmin, fmax]
    %
    % Output:
    %   SPKF: Spike-field coherence for each neuron-channel pair [NeuronID x Channel x Frequency]

    [nTrials, nChans, nTime] = size(LFP);
    [nTrials2, nNeurons, nTime2] = size(SPK);

    % Ensure dimensions match
    assert(nTrials == nTrials2, 'Number of trials in LFP and SPK must match');
    assert(nTime == nTime2, 'Time points must match in LFP and SPK');

    % Define frequency range
    fmin = fRange(1);
    fmax = fRange(2);
    
    % Frequency vector
    freqs = fmin:(fmax - fmin) / 100:fmax;

    % Initialize output coherence matrix [Neuron x Channel x Frequency]
    SPKF = zeros(nNeurons, nChans, length(freqs));

    % Loop over each neuron and each channel
    for neuron = 1:nNeurons
        for chan = 1:nChans
            % Extract LFP and Spike for current channel and neuron
            lfpData = squeeze(LFP(:, chan, :));  % [Trial x Time]
            spkData = squeeze(SPK(:, neuron, :)); % [Trial x Time]

            % Initialize accumulators for spectra
            Sxy_sum = 0;
            Sxx_sum = 0;
            Syy_sum = 0;

            % Loop through each trial
            for trial = 1:nTrials
                % Extract single trial LFP and spike data
                lfp_trial = lfpData(trial, :);
                spk_trial = spkData(trial, :);

                % Compute the cross-spectrum between LFP and spikes
                [Sxy, f] = cpsd(lfp_trial, spk_trial, [], [], freqs, Fs);

                % Compute the power spectra
                [Sxx, ~] = pwelch(lfp_trial, [], [], freqs, Fs);
                [Syy, ~] = pwelch(spk_trial, [], [], freqs, Fs);

                % Accumulate the spectra over trials
                Sxy_sum = Sxy_sum + Sxy;
                Sxx_sum = Sxx_sum + Sxx;
                Syy_sum = Syy_sum + Syy;
            end

            % Average over trials
            Sxy_mean = Sxy_sum / nTrials;
            Sxx_mean = Sxx_sum / nTrials;
            Syy_mean = Syy_sum / nTrials;

            % Compute coherence
            coherence = abs(Sxy_mean).^2 ./ (Sxx_mean .* Syy_mean);
            SPKF(neuron, chan, :) = coherence;
        end
    end
end
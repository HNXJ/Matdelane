function output = j4DTimeFrequencyResponse(input)

    [sessions, trials, channels, ~] = size(input);
    
    % Set parameters
    freqLimits = [0 150];
    timeResolution = 0.4;
    overlap = 97;
    
    % Preallocate output matrix
    [~, f, t] = pspectrum(squeeze(input(1,1,1,:)), 'FrequencyLimits', freqLimits, ...
                          'TimeResolution', timeResolution, 'OverlapPercent', overlap);
    output = zeros(sessions, trials, channels, length(f), length(t));
    
    % Use parfor for parallel processing
    parfor s = 1:sessions
        for tr = 1:trials
            for ch = 1:channels
                [tfr, ~, ~] = pspectrum(squeeze(input(s,tr,ch,:)), 'FrequencyLimits', freqLimits, ...
                                        'TimeResolution', timeResolution, 'OverlapPercent', overlap);
                output(s,tr,ch,:,:) = tfr;
            end
        end
    end

end
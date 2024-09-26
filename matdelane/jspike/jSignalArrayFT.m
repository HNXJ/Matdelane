function ysignal = jSignalArrayFT(xsignal, channelTag)

    Fs = 1000; % Sampling rate in Hz
    Ls = size(xsignal, 3); % Length of signal
    
    nChannels = size(xsignal, 2);
    nTrials = size(xsignal, 1);
    
    for iChannel = 1:nChannels
    
        tempsig = squeeze(xsignal(:, iChannel, :));
    
        for iTrial = 1:nTrials
    
            ysignal.trial{iTrial}(iChannel, :) = tempsig(iTrial, :); 
            ysignal.time{iTrial} = (1:Ls)/Fs;
            ysignal.label{iChannel} = [channelTag '-' num2str(iChannel)]; 
            ysignal.fsample = Fs; 
    
        end
    
    end
    
    ysignal = ft_checkdata(ysignal, 'datatype', 'raw', 'feedback', 'yes');

end
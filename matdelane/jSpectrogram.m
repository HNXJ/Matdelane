function [y, t, f] = jSpectrogram(X, fs, fmax, freqW, timeW, overlap, fcnt, tcnt, smoothKernel, plotFlag, logFlag)

%{

[y, t, f] = jSpectrogram(X, fs, fmax, fcnt, tcnt, freqW, timeW, overlap, smoothKernel, plotFlag)

>Input args 
x : Signal (<1 x N> numeric array)
fs : sampling rate (scalar, Hz)
fmax : Frequency upper bound (scalar, Hz)
fcnt : Frequency bin count (numeric scalar)
tcnt : Time bin count (scalar, ms)
freqW : Frequency bin width for the time-frequency response (scalar, Hz)
timeW : Time bin width for the time-frequency response (scalar, ms)
overlap : Overlap between time windows (scalar, ms)
smoothKernel : Gaussian smoothing kernel window size (scalar, ms)
plotFlag : Plotting status (boolean)

>Output args
y : Time-frequency power amplitudes (<tcnt x fcnt> numeric array)
t : Time bins (<1 x tcnt> numeric array)
f : Frequency bins (<1 x fcnt> numeric array) 

%}


    m = size(X, 1);
    xfactor = fs / 1000;
    tmax = size(X, 2) * xfactor;

    if ~exist('fs', 'var')

        fs = 1000;

    end

    if ~exist('fmax', 'var')

        fmax = 100;

    end

    if ~exist('freqW', 'var')

        freqW = 1;

    end

    if ~exist('timeW', 'var')

        timeW = 100;

    end

    if ~exist('overlap', 'var')

        overlap = 90;

    end

    if ~exist('fcnt', 'var')

        fcnt = floor(fmax / freqW);

    end

    tW = ((timeW - overlap)*xfactor);

    if ~exist('tcnt', 'var')

        tcnt = floor(tmax / tW);

    end

    if ~exist('smoothKernel', 'var')

        smoothKernel = 10;

    end

    if ~exist('plotFlag', 'var')

        plotFlag = 0;

    end

    if ~exist('logFlag', 'var')

        logFlag = 1;

    end

    tWx = (timeW)*(fs/1000);
    tB = tcnt;
    fB = fcnt;

    kernelSize = ceil(xfactor*smoothKernel);
    y = zeros(m, tB, fB);

    parfor i = 1:m

        for j = 1:tB

            lt = max(j*tW - tWx, 1);
            rt = max(j*tW, 1);
            tK = lt:rt;
            tempX = X(i, tK);
            tempT = exp(-abs(linspace(-1.0, 1.0, kernelSize).^2));
            tempX = conv(tempX, tempT/sum(tempT), "same");
            tempF = jSpectrum(tempX, fs, fmax, fB, 0);
            y(i, j, :) = smooth(tempF, 10);

        end

    end

    t = linspace(0, tmax, tB);
    f = linspace(0, fmax, fB);
    sG = squeeze(mean(y, 1));
    y = sG';
    y(:, 1:50) = y(:, mod(0:14, 10)+6);
    yx = mean(y(:, 1:20), 2);

    for i = 1:size(y, 2)

        y(:, i) = y(:, i) ./ yx;

    end


    % y(:, 1:25) = 1;

    if plotFlag

        % if logFlag
        % 
        %     y = log(y);
        % 
        % end

        figure('Position', [0, 0, 1700, 1400]);
        subplot(1, 1, 1);
        imagesc(y-1, "XData", t, "YData", f);
        clim([-0.4 0.4]);
        
        xlabel("Time (ms)");
        ylabel("Freq (Hz)");
        % colormap("jet");
        colorbar();
        set(gca,'YDir','normal');

        sgtitle("Spectrogram");

    end

end
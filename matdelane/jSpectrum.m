function [y, f] = jSpectrum(x, fs, fmax, fcnt, smoothing)

%{
[y, f] = jSpectrum(x, fs, fmax, fcnt, smoothing)

>Input args 
x : Signal (<1 x N> numeric array)
fs : sampling rate (scalar, Hz)
fmax : Frequency upper bound (scalar, Hz)
fcnt : Frequency bin count (numeric scalar)
smoothing : Gaussian low-pass filter (numeric scalar)

>Output args
y : Power amplitudes (<1 x fcnt> numeric array)
f : Frequency bins (<1 x fcnt> numeric array) 

%}

    N = length(x);

    if ~exist('fs', 'var')

        fs = 1000;
        disp("No sampling rate was specificed. Default Fs = 1000Hz ");

    end

    if ~exist('fmax', 'var')

        fmax = fs/2;

    end

    if ~exist('fcnt', 'var')

        fcnt = N;

    end

    if ~exist('smoothing', 'var')

        smoothing = 0;

    end

    if fmax*2 > fs

        fmax = fs/2;

    end

    t1 = linspace(0, fs/2, ceil(N/2)-1);
    t2 = linspace(0, fmax, fcnt);

    if smoothing

        x = smooth(x, smoothing);

    end

    p1 = fft(x);
    p2 = abs(p1(2:ceil(N/2)));
    y = interp1(t1, p2, t2);
    f = t2;

end
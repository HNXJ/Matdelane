function [phi, frq] = jInstantSpectrum(x, Fs, phaseSmooth)

    % x : 1-D signal
    % Fs : Sampling rate
    % phaseSmooth (opt) : width of phase smoothing kernel

    if ~exist("phaseSmooth", "var")

        phaseSmooth = 20;

    end

    x = smooth(x, 10);

    N = length(x);
    dT = 2*pi / Fs;

    hx = hilbert(x);
    sa = x + 1i*hx;
    phi = phase(sa);
    % phi = smooth(phi, phaseSmooth);
    frq = diff(phi)/(dT);
    frq(1) = frq(2);
    frq(end) = frq(end-1);

end
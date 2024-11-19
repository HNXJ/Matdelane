function [phi, frq] = jInstantSpectrum(x, Fs, psmooth)

    % x : 1-D signal
    % Fs : Sampling rate
    % phaseSmooth (opt) : width of phase smoothing kernel

    if ~exist("psmooth", "var")

        psmooth = 20;

    end

    x = smooth(x, psmooth);
    dT = 2*pi / Fs;

    hx = hilbert(x);
    sa = x + 1i*hx;
    phi = phase(sa);
    frq = diff(phi)/(dT);
    frq(1) = frq(2);
    frq(end) = frq(end-1);

end
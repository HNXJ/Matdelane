function [pgx, xinfo] = jCalculateTFRbyChannel(xlfp, channel_in_layer)

    freqlims = [0 200];
    xinfo.freqlims = freqlims;
    leakage = 0.85;
    xinfo.leakage = leakage;

    freqres = 5.0;
    xinfo.freqres = freqres;
    overlap = 95;
    xinfo.overlap = overlap;
    
    channel_in_layer_selected = channel_in_layer;
    y = squeeze(mean(mean(xlfp{1}, 1), 2));tpre = 500;
    [p1, xinfo.fmap, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", xinfo.freqlims, "OverlapPercent", xinfo.overlap, "FrequencyResolution", xinfo.freqres, "Leakage", xinfo.leakage);
    xinfo.tmap = (t1 - tpre/1000)*1000;% TFR Kaiser's time window offset shift = t(1)*2
    
    pgx = cell(12, 4);
    
    for ik = 1:12
        
        xG1d = xlfp{ik}(:, channel_in_layer_selected.deep, :);
        xG1m = xlfp{ik}(:, channel_in_layer_selected.mid, :);
        xG1s = xlfp{ik}(:, channel_in_layer_selected.sup, :);
        xG1f = xlfp{ik}(:, channel_in_layer_selected.goodch, :);
        
        TR = size(xG1d, 1);
    
        pgxcd = zeros([TR, size(p1)]);
        pgxcm = zeros([TR, size(p1)]);
        pgxcs = zeros([TR, size(p1)]);
        pgxcf = zeros([TR, size(p1)]);
        
        parfor jk = 1:TR % trials
        
            dCH = size(xG1d, 2);

            for kk = 1:dCH

                [p1temp, ~, ~] = pspectrum(squeeze(xG1d(jk, kk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
                pgxcd(jk, :, :) = p1temp + squeeze(pgxcd(jk, :, :));

            end

            mCH = size(xG1m, 2);
            
            for kk = 1:mCH

                [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, kk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
                pgxcm(jk, :, :) = p1temp + squeeze(pgxcm(jk, :, :));

            end
              
            sCH = size(xG1s, 2);
            
            for kk = 1:sCH

                [p1temp, ~, ~] = pspectrum(squeeze(xG1s(jk, kk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
                pgxcs(jk, :, :) = p1temp + squeeze(pgxcs(jk, :, :));

            end      
                    
            fCH = size(xG1f, 2);
            
            for kk = 1:fCH

                [p1temp, ~, ~] = pspectrum(squeeze(xG1f(jk, kk, :)), 1000, "spectrogram", "FrequencyLimits", freqlims, "FrequencyResolution", freqres, "OverlapPercent", overlap, "Leakage", leakage);
                pgxcf(jk, :, :) = p1temp + squeeze(pgxcf(jk, :, :));

            end
    
            if mod(jk, 20) == 0
    
                fprintf(num2str(jk));
    
            end
        
        end
    
        pgx{ik, 1} = pgxcd;
        pgx{ik, 2} = pgxcm;
        pgx{ik, 3} = pgxcs;
        pgx{ik, 4} = pgxcf;
    
        disp(" >cond : " + num2str(ik));
    
    end

    xinfo.tbands = cell(1, 5);
    xinfo.tbands{1} = find(xinfo.tmap > -250, 1):find(xinfo.tmap > -50, 1);
    xinfo.tbands{2} = find(xinfo.tmap > 0, 1):find(xinfo.tmap > 1000, 1);
    xinfo.tbands{3} = find(xinfo.tmap > 1031, 1):find(xinfo.tmap > 2031, 1);
    xinfo.tbands{4} = find(xinfo.tmap > 2062, 1):find(xinfo.tmap > 3062, 1);
    xinfo.tbands{5} = find(xinfo.tmap > 3093, 1):find(xinfo.tmap > 4093, 1);
    
    nt_temp = length(xinfo.tbands{3});
    
    for ik = 2:5
    
        xinfo.tbands{ik} = xinfo.tbands{ik}(1:nt_temp);
    
    end

    
    xinfo.fbands = cell(1, 5);
    
    xinfo.fbands{1} = find(xinfo.fmap > 2, 1):find(xinfo.fmap > 7, 1);
    xinfo.fbands{2} = find(xinfo.fmap > 8, 1):find(xinfo.fmap > 12, 1);
    xinfo.fbands{3} = find(xinfo.fmap > 14, 1):find(xinfo.fmap > 30, 1);
    xinfo.fbands{4} = find(xinfo.fmap > 32, 1):find(xinfo.fmap > 80, 1);
    xinfo.fbands{5} = find(xinfo.fmap > 80, 1):find(xinfo.fmap >= max(xinfo.fmap), 1);

    return;

end
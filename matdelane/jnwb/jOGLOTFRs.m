function [pgx, pgxcnt] = jOGLOTFRs(xG)
    
    y = squeeze(mean(mean(xG{1, 1}{1}, 1), 2));
    frqlims = [0 250];
    timeres = 0.2;
    overlap = 95;
    
    [p1, ~, ~] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", frqlims, "TimeResolution", timeres, "OverlapPercent", overlap);
    sessionsN = size(xG{1, 1}, 2);
    pgx = zeros([sessionsN, 12, 128, size(p1)]);
    pgxcnt = zeros([sessionsN, 12, 128, 2]);
    
    for ik = 1:12 % conditions
    
        xG1 = xG{1, ik};
        xGc1 = xG{2, ik};
    
        for jk = 1:size(xG1, 2) % sessions
    
            xG2 = xG1{jk};
            xGc2 = xGc1{jk}.rfdata;
    
            qRF = xGc2(:, 4) > 1.0;
            xRF = xGc2(:, 1);
            yRF = xGc2(:, 2);
            rRF = sqrt(xRF.^2 + yRF.^2);
    
            fovealch = rRF < 2.0;
            pfovealch = (rRF >= 2.0) & (rRF <= 5.0);
            nfovealch = rRF > 5.0;
    
            for kk = 1:size(xG2, 1) % trial
    
                for lk = 1:size(xG2, 2) % channel
    
                    pgxcnt(jk, ik, lk, 1) = pgxcnt(jk, ik, lk, 1) + 1;
    
                    if fovealch(lk)
    
                        pgxcnt(jk, ik, lk, 2) = 1;
    
                    elseif pfovealch(lk)
    
                        pgxcnt(jk, ik, lk, 2) = 2;
    
                    elseif nfovealch(lk)
    
                        pgxcnt(jk, ik, lk, 2) = 3;
    
                    else
    
                        pgxcnt(jk, ik, lk, 2) = -1;
    
                    end
                
                end
    
            end
    
        end
    
        disp("cond : " + num2str(ik));
    
    end
    
    for ik = 1:12 % conditions
    
        xG1 = xG{1, ik};
    
        for jk = 1:size(xG1, 2) % sessions
    
            xG2 = xG1{jk};
    
            parfor lk = 1:size(xG2, 2) % channel
    
                y = squeeze(mean(xG2(:, lk, :), 1, "omitnan"));
                
                if sum(y) > 0
                    [pgx(jk, ik, lk, :, :), ~, ~] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", frqlims, "TimeResolution", timeres, "OverlapPercent", overlap);
                end
            
            end
    
            fprintf("ses" + num2str(jk));
    
        end
    
        disp("cond : " + num2str(ik));
    
    end
    
end
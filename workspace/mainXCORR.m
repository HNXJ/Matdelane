%% Stim

fbandlabels = {"Theta[2.5-7Hz]", "Alpha[8-12Hz]", "Beta[13-30Hz]", "L-Gamma[32-80Hz]", "H-Gamma[80-200Hz]"};
areanames = {"V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"};

figure("Position", [0 0 1200 2200]);

for fbandx = 1:4

    imx1 = corr(gcorr_stim{fbandx}', gcorr_stim{5}', "Type", "Spearman");

    gcorr_a = zeros(11, 11);

    for ik = 1:Nfiles

        if ik == 1

            lcntx = 0;

        else

            lcntx = areacnt(ik-1);

        end
        
        for jk = 1:Nfiles

            if jk == 1
    
                rcntx = 0;
    
            else
    
                rcntx = areacnt(jk-1);
    
            end

            gcorr_a(ik, jk) = mean(imx1(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");

        end

    end

    % imx1(abs(imx1) < 0.5);

    subplot(2, 2, fbandx);

    gcorr_a(abs(gcorr_a) < 0.01) = 0;
    gcorr_a = gcorr_a*100;

    imagesc(gcorr_a);

    for ik = 1:Nfiles

        for jk = 1:Nfiles
            
            if gcorr_a(ik, jk)
                
                text(jk-.3, ik, num2str(gcorr_a(ik, jk), 3));

            else

                text(jk-.3, ik, "n.s");

            end
            
        end

    end
    clim([-100 100]);
    
    title("Spectral power corr" + fbandlabels{fbandx} + " to hGamma");
    xticks(1:11);
    xticklabels(areanames);
    yticks(1:11);
    yticklabels(areanames);
    colormap("winter");

end


sgtitle("stim time xCorr (Sxx/Sxy)"); % TODO make |r| < .25 0
fname = "allareatfr_hcorr" + "xstim";
print(gcf,'-vector','-dsvg', fname +".svg");

%% Omission

fbandlabels = {"Theta[2.5-7Hz]", "Alpha[8-12Hz]", "Beta[13-30Hz]", "L-Gamma[32-80Hz]", "H-Gamma[80-200Hz]"};
areanames = {"V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"};

figure("Position", [0 0 1200 2200]);

for fbandx = 1:4

    imx1 = corr(gcorr_omxn{fbandx}', gcorr_omxn{5}', "Type", "Spearman");

    gcorr_a = zeros(11, 11);

    for ik = 1:Nfiles

        if ik == 1

            lcntx = 0;

        else

            lcntx = areacnt(ik-1);

        end
        
        for jk = 1:Nfiles

            if jk == 1
    
                rcntx = 0;
    
            else
    
                rcntx = areacnt(jk-1);
    
            end

            gcorr_a(ik, jk) = mean(imx1(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");

        end

    end

    % imx1(abs(imx1) < 0.5);

    subplot(2, 2, fbandx);

    gcorr_a(abs(gcorr_a) < 0.01) = 0;
    gcorr_a = gcorr_a*100;

    imagesc(gcorr_a);

    for ik = 1:Nfiles

        for jk = 1:Nfiles
            
            if gcorr_a(ik, jk)
                
                text(jk-.3, ik, num2str(gcorr_a(ik, jk), 3));

            else

                text(jk-.3, ik, "n.s");

            end
            
        end

    end
    clim([-50 50]);
    
    title("Spectral power corr" + fbandlabels{fbandx} + " to hGamma");
    xticks(1:11);
    xticklabels(areanames);
    yticks(1:11);
    yticklabels(areanames);
    % colormap("winter");

end


sgtitle("Ox time xCorr (Sxx/Sxy)"); % TODO make |r| < .25 0
fname = "allareatfr_hcorr" + "xomxn";
print(gcf,'-vector','-dsvg', fname +".svg");

%%

fbandlabels = {"Theta[2.5-7Hz]", "Alpha[8-12Hz]", "Beta[13-30Hz]", "L-Gamma[32-80Hz]", "H-Gamma[80-200Hz]"};
areanames = {"V1", "V2", "V3d", "V3a", "V4", "MT", "MST", "TEO", "FST", "FEF", "PFC"};

figure("Position", [0 0 1200 2200]);

for fbandx = 1:4

    imx1 = corr(gcorr_stim{fbandx}', gcorr_stim{5}', "Type", "Spearman");
    imx2 = corr(gcorr_omxn{fbandx}', gcorr_omxn{5}', "Type", "Spearman");

    gcorr_a = zeros(11, 11);
    gcorr_b = zeros(11, 11);

    for ik = 1:Nfiles

        if ik == 1

            lcntx = 0;

        else

            lcntx = areacnt(ik-1);

        end
        
        for jk = 1:Nfiles

            if jk == 1
    
                rcntx = 0;
    
            else
    
                rcntx = areacnt(jk-1);
    
            end

            gcorr_a(ik, jk) = mean(imx1(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");
            gcorr_b(ik, jk) = mean(imx2(lcntx+1:areacnt(ik), rcntx+1:areacnt(jk)), "all");

        end

    end

    % imx1(abs(imx1) < 0.5);

    subplot(2, 2, fbandx);

    gcorr_c = 100*(-gcorr_b.^2 + gcorr_a.^2);
    gcorr_c(isnan(gcorr_c)) = 0;
    gcorr_c(abs(gcorr_c) < 1) = 0;

    imagesc(gcorr_c);

    for ik = 1:Nfiles

        for jk = 1:Nfiles

            if gcorr_c(ik, jk)
                
                text(jk-.3, ik, num2str(gcorr_c(ik, jk), 3));

            else

                text(jk-.3, ik, "n.s");

            end
        end

    end

    clim([-20 20]);
    
    title("Spectral power corr ratio  Stim - Oxm" + fbandlabels{fbandx} + " to hGamma");
    xticks(1:11);
    xticklabels(areanames);
    yticks(1:11);
    yticklabels(areanames);
    colormap("winter");

end

sgtitle("Corr change ( stim - omission )"); % TODO make |r| < .25 0
fname = "allareatfr_hcorr" + "_r2change";
print(gcf,'-vector','-dsvg', fname +".svg");

%%
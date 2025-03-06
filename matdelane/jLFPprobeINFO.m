function jLFPprobeINFO(x)

    imglfp1 = squeeze(mean(x, 1));

    figure;
    subplot(2, 1, 1);
    imagesc(imglfp1);
    xlabel("Time (ms)");
    ylabel("Channel");
    title("LFP average by trial");
    
    corrx = corr(imglfp1', "Type", "Spearman");
    subplot(2, 2, 3);
    imagesc(corrx);
    xlabel("Channel");
    ylabel("Channel");
    title("Corr of LFP average by trial");
    colorbar;
    
    subplot(2, 2, 4);
    imxcorr = mean(corr(corrx), 1);
    stem(imxcorr);
    xlabel("Channel");
    ylabel("Spearman deviation index");
    title("Channel grand corr. (corr of corr)");
    
end
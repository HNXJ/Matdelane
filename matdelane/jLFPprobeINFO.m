function jLFPprobeINFO(x, chx)

    imglfp1 = squeeze(mean(x, 1));

    figure;
    subplot(4, 1, 1);
    imagesc(imglfp1, "YData", chx);
    xlabel("Time (ms)");
    ylabel("Channel");
    title("LFP average by trial");

    imglfp2 = smoothdata2(imglfp1, "gaussian", {10, 20});
    imgcsd1 = (-imglfp2(1:end-2, :) + 2*imglfp2(2:end-1, :) - imglfp2(3:end, :));
    imgcsd1 =  smoothdata2(imgcsd1, "gaussian", {7, 20});
    % imgcsd1 = (imgcsd1);

    subplot(4, 1, 2);
    imagesc(imgcsd1, "YData", chx);
    xlabel("Time (ms)");
    ylabel("Channel");%clim([-1 1]);
    title("CSD average by trial");

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
    hold("on");

    % rectangle("Position", [1 2 3 4], "FaceColor", "b", "FaceAlpha", 0.25);
    % rectangle("Position", [1 2 3 4], "FaceColor", "b", "FaceAlpha", 0.25);
    % rectangle("Position", [1 2 3 4], "FaceColor", "b", "FaceAlpha", 0.25);
    
    xlabel("Channel");
    ylabel("Spearman deviation index");
    title("Channel grand corr. (corr of corr)");
    
end
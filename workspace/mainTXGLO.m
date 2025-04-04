%% Concatenate OGLO analysis


q1 = load("OGLOobj\sub-C31o_ses-230816PFC.mat", "obj").obj;
q2 = load("OGLOobj\sub-C31o_ses-230818PFC.mat", "obj").obj;

%%


q1.jTFRplot(9, 4, q1.tbands{2}(end-10:end-5));
q2.jTFRplot(9, 4, q2.tbands{2}(end-10:end-5));

%% TFR barplot



%%

tcond1 = 9;
layerid = 4;
tbaseline = q1.tbands{2}(end-10:end-5);
txlims = [obj.tmap(1) obj.tmap(end)];

figure;
subplot(2, 1, 1);
tfr1 = squeeze(mean(q1.pgx{tcond1, layerid}, 1)) + squeeze(mean(q2.pgx{tcond1, layerid}, 1));
% tfr1e = squeeze(std(q1.pgx{tcond1, layerid}));

for ik = 1:size(tfr1, 1)

    tfr1(ik, :) = tfr1(ik, :) / mean(tfr1(ik, tbaseline));

end

imagesc(10*log10(tfr1), "XData", q1.tmap, "YData", q1.fmap);
% ylim([0 20])
xlim(txlims);
set(gca, "YDir", "normal");
hold("on");
xline(0);
xline(1031);
xline(2062);
xline(3093);
xlabel("Time (ms)");
ylabel("Frequency (Hz)");
title("Power change from baseline (dB)");
colorbar;

for fband = 1:5

    yline(q1.fmap(q1.fbands{fband}(1)), "Color", [1 0 0]);

end

xlabel("Time (ms)");
subplot(2, 1, 2);
cls = zeros(5, 3);
cls(1, :) = [0 0 1];
cls(2, :) = [1 0 0];
cls(3, :) = [1 0.5 0];
cls(4, :) = [1 0 1];
cls(5, :) = [0 0.8 0.4];

for fband = 1:5

    tfr1 = squeeze(mean(q1.pgx{tcond1, layerid}(:, q1.fbands{fband}, :), 1)) + squeeze(mean(q2.pgx{tcond1, layerid}(:, q2.fbands{fband}, :), 1));
    tfr1 = squeeze(mean(tfr1, 1));
    tfr1 = smooth(tfr1, 10);
    baselinecrx = mean(tfr1(tbaseline));

    tfr1 = tfr1 / baselinecrx;

    ne1 = length(q1.fbands{fband}); % spectral res. bands
    ne2 = size(q1.pgx{tcond1, layerid}, 1); % trial count
    ne3 = size(q1.x{tcond1}, 2); % channel count

    tfr1e = squeeze(std(q1.pgx{tcond1, layerid}(:, q1.fbands{fband}, :) / baselinecrx));
    tfr1e = squeeze(mean(tfr1e, 1)) / sqrt(ne1 + ne2 + ne3);
    tfr1e = smooth(tfr1e, 10);

    cl = cls(fband, :);
    y1s = 10*log10(tfr1);
    plot(q1.tmap, y1s, "DisplayName", q1.fbandlabels(fband), "LineWidth", 0.5, "Color", cl);
    hold("on");

    stx = tfr1 + 2*tfr1e;
    sty = tfr1 - 2*tfr1e;
    stx(stx <= 0) = 1e-2;
    sty(sty <= 0) = 1e-2;

    stx = 10*log10(stx);
    sty = 10*log10(sty);

    patch([q1.tmap; q1.tmap(end:-1:1)], [stx; sty(end:-1:1)], cl, "EdgeColor", "none", "FaceColor", cl, "FaceAlpha", 0.5, "HandleVisibility", "off");

    % ylim([0 20])
    xlim(txlims);

    xline(0, HandleVisibility="off");
    xline(1031, HandleVisibility="off");
    xline(2062, HandleVisibility="off");
    xline(3093, HandleVisibility="off");
    xlabel("Times (ms)");
    ylabel("Power change (dB)");
    title("Power change from fixation baseline (+-2Se)");

end

colorbar;
sgtitle(q1.areainf + q1.condinflabel(tcond1) + "/" + q1.layeridlabel(layerid));
legend;
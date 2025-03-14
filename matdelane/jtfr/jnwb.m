classdef jnwb < handle

    % JNWB: Time-Frequency analysis on signals; Passive omission task
    %  LFP, SPK (SUA and MUA)
    
    properties

        nwbFile = "_file.nwb";
        areainf = "_area_/";
        layeridlabel = ["deep", "mid", "sup", "all"];
        condinflabel = ["AAAB", "AXAB", "AAXB", "AAAX", "BBBA", "BXBA", "BBXA",...
                "BBBX", "RRRR", "RXRR", "RRXR", "RRRX"];

        c = {};
        x = {};
        cm = {};
        xm = {};
        cs = {};
        xs = {};
        
        tmapsignal = 0;
        tpre = 0;
        tpost = 0;
        probeid = 0;

        freqres = 0;
        convflag = 0;
        tbands = 0;
        fbands = 0;

        pgx = {};
        leakage = 0;
        freqlims = 0;
        overlap = 0;

        fmap = 0;
        tmap = 0;
        tbandlabels = ["Fix", "S1d1", "S2d2", "S3d3", "S4d4"];
        fbandlabels = ["Theta(2-7Hz)", "Alpha(8-12Hz)", "Beta(14-30Hz)", "GammaL(32-80Hz)", "GammaH(80+Hz)"];

    end
    
    methods

        function obj = jnwb(nwbFile, areainf, tpre, tpost, probeid, convFlag)

            disp(">Loading : " + nwbFile);

            obj.areainf = areainf;
            obj.nwbFile = nwbFile;
            obj.tpre = tpre;
            obj.tpost = tpost;

            obj.probeid = probeid;
            obj.convflag = convFlag;
            nwb = nwbRead(nwbFile);

            [obj.c, obj.x] = jOGLOSignals(nwb, "omission_glo_passive", tpre, tpost, probeid, "lfp");
            [obj.cm, obj.xm] = jOGLOSignals(nwb, "omission_glo_passive", tpre, tpost, probeid, "muae");
            [obj.cs, obj.xs] = jOGLOUnits(nwb, "omission_glo_passive", tpre, tpost, convFlag);

            obj.tmapsignal = linspace(-obj.tpre, obj.tpost, size(obj.x{1}, 3));

        end

        function obj = tmapsignalset(obj)

            obj.tmapsignal = linspace(-obj.tpre, obj.tpost, size(obj.x{1}, 3));

        end

        function jCalcTFR(obj, channel_in_layer)

            obj.freqlims = [0 200];
            obj.leakage = 0.85;
            obj.freqres = 5.0;
            obj.overlap = 95;
            
            channel_in_layer_selected = channel_in_layer;
            y = squeeze(mean(mean(obj.x{1}, 1), 2));
            [p1, obj.fmap, t1] = pspectrum(y, 1000, "spectrogram", "FrequencyLimits", obj.freqlims, "OverlapPercent", obj.overlap, "FrequencyResolution", obj.freqres, "Leakage", obj.leakage);
            obj.tmap = (t1 - obj.tpre/1000)*1000;% TFR Kaiser's time window offset shift = t(1)*2
            
            obj.pgx = cell(12, 4);
            
            for ik = 1:12
                
                xG1d = squeeze(mean(obj.x{ik}(:, channel_in_layer_selected.deep, :), 2));
                xG1m = squeeze(mean(obj.x{ik}(:, channel_in_layer_selected.mid, :), 2));
                xG1s = squeeze(mean(obj.x{ik}(:, channel_in_layer_selected.sup, :), 2));
                xG1f = squeeze(mean(obj.x{ik}(:, channel_in_layer_selected.goodch, :), 2));
                
                TR = size(xG1d, 1);
            
                pgxcd = zeros([TR, size(p1)]);
                pgxcm = zeros([TR, size(p1)]);
                pgxcs = zeros([TR, size(p1)]);
                pgxcf = zeros([TR, size(p1)]);
                
                parfor jk = 1:TR % trials
                
                    [p1temp, ~, ~] = pspectrum(squeeze(xG1d(jk, :)), 1000, "spectrogram", "FrequencyLimits", obj.freqlims, "FrequencyResolution", obj.freqres, "OverlapPercent", obj.overlap, "Leakage", obj.leakage);
                    pgxcd(jk, :, :) = p1temp;
                
                    [p1temp, ~, ~] = pspectrum(squeeze(xG1m(jk, :)), 1000, "spectrogram", "FrequencyLimits", obj.freqlims, "FrequencyResolution", obj.freqres, "OverlapPercent", obj.overlap, "Leakage", obj.leakage);
                    pgxcm(jk, :, :) = p1temp;
                
                    [p1temp, ~, ~] = pspectrum(squeeze(xG1s(jk, :)), 1000, "spectrogram", "FrequencyLimits", obj.freqlims, "FrequencyResolution", obj.freqres, "OverlapPercent", obj.overlap, "Leakage", obj.leakage);
                    pgxcs(jk, :, :) = p1temp;
                     
                    [p1temp, ~, ~] = pspectrum(squeeze(xG1f(jk, :)), 1000, "spectrogram", "FrequencyLimits", obj.freqlims, "FrequencyResolution", obj.freqres, "OverlapPercent", obj.overlap, "Leakage", obj.leakage);
                    pgxcf(jk, :, :) = p1temp;
            
                    if mod(jk, 20) == 0
            
                        fprintf(num2str(jk));
            
                    end
                
                end
            
                obj.pgx{ik, 1} = pgxcd;
                obj.pgx{ik, 2} = pgxcm;
                obj.pgx{ik, 3} = pgxcs;
                obj.pgx{ik, 4} = pgxcf;
            
                disp(" >cond : " + num2str(ik));
            
            end

            obj.tbands = cell(1, 5);
            obj.tbands{1} = find(obj.tmap > -250, 1):find(obj.tmap > -50, 1);
            obj.tbands{2} = find(obj.tmap > 0, 1):find(obj.tmap > 1000, 1);
            obj.tbands{3} = find(obj.tmap > 1031, 1):find(obj.tmap > 2031, 1);
            obj.tbands{4} = find(obj.tmap > 2062, 1):find(obj.tmap > 3062, 1);
            obj.tbands{5} = find(obj.tmap > 3093, 1):find(obj.tmap > 4093, 1);
            
            nt_temp = length(obj.tbands{3});
            
            for ik = 2:5
            
                obj.tbands{ik} = obj.tbands{ik}(1:nt_temp);
            
            end

            
            obj.fbands = cell(1, 5);
            
            obj.fbands{1} = find(obj.fmap > 2, 1):find(obj.fmap > 7, 1);
            obj.fbands{2} = find(obj.fmap > 8, 1):find(obj.fmap > 12, 1);
            obj.fbands{3} = find(obj.fmap > 14, 1):find(obj.fmap > 30, 1);
            obj.fbands{4} = find(obj.fmap > 32, 1):find(obj.fmap > 80, 1);
            obj.fbands{5} = find(obj.fmap > 80, 1):find(obj.fmap >= max(obj.fmap), 1);

            return;

        end

        function jMUAplot(obj, condid, timewindow)

            if ~exist("timewindow", "var")

                obj.tmapsignalset();
                timewindow = [obj.tmapsignal(1), obj.tmapsignal(end)];

            end

            figure;
            
            imxm = squeeze(mean(obj.xm{condid}, 1));
            imxm = squeeze(mean(imxm, 1));
            imxm = (imxm - mean(imxm)) / std(imxm);
            plot(linspace(-500, 4250, 4750), imxm, "DisplayName", obj.areainf);
            
            hold("on");
            xline(0, HandleVisibility="off");
            xline(1031, HandleVisibility="off");
            xline(2062, HandleVisibility="off");
            xline(3093, HandleVisibility="off");
            
            title("MUAenv/Zsc/" + obj.condinflabel(condid));
            xlabel("Time (ms)");
            ylabel("Z-score");
            xlim(timewindow);
            
            legend;

        end

        function jSUAplot(obj, condid, timewindow, units)

            if ~exist("timewindow", "var")

                obj.tmapsignalset();
                timewindow = [obj.tmapsignal(1), obj.tmapsignal(end)];

            end

            if ~exist("units", "var")

                units = 1:size(obj.xs{condid}, 2);

            end

            figure;
            
            imxm = squeeze(mean(obj.xs{condid}(:, units, :), 1));
            imxm = squeeze(mean(imxm, 1));
            imxm = (imxm - mean(imxm)) / std(imxm);
            plot(linspace(-500, 4250, 4750), imxm, "DisplayName", obj.areainf);
            
            hold("on");
            xline(0, HandleVisibility="off");
            xline(1031, HandleVisibility="off");
            xline(2062, HandleVisibility="off");
            xline(3093, HandleVisibility="off");
            
            title("SUAenv/Zsc/" + obj.condinflabel(condid));
            xlabel("Time (ms)");
            ylabel("Z-score");
            xlim(timewindow);
            
            legend;

        end

        function jTFRplot(obj, tcond1, layerid, tbaseline, txlims)

            if ~exist("tbaseline", "var")
            
                tbaseline = obj.tbands{1};

            end

            if ~exist("txlims", "var")
            
                txlims = [obj.tmap(1) obj.tmap(end)];

            end
            
            figure;
            subplot(2, 1, 1);
            tfr1 = squeeze(mean(obj.pgx{tcond1, layerid}, 1));
            
            for ik = 1:size(tfr1, 1)
            
                tfr1(ik, :) = tfr1(ik, :) / mean(tfr1(ik, tbaseline));
            
            end
            
            imagesc(10*log10(tfr1), "XData", obj.tmap, "YData", obj.fmap);
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
            
                yline(obj.fmap(obj.fbands{fband}(1)), "Color", [1 0 0]);
            
            end
            
            xlabel("Time (ms)");
            subplot(2, 1, 2);
            
            for fband = 1:5
            
                tfr1 = squeeze(mean(obj.pgx{tcond1, layerid}(:, obj.fbands{fband}, :), 1));
                tfr1 = squeeze(mean(tfr1, 1));
                tfr1 = smooth(tfr1, 20);
                tfr1 = tfr1 / mean(tfr1(tbaseline));
                plot(obj.tmap, 10*log10(tfr1), "DisplayName", obj.fbandlabels(fband), "LineWidth", 2);
                % ylim([0 20])
                xlim(txlims);
                hold("on");
                xline(0, HandleVisibility="off");
                xline(1031, HandleVisibility="off");
                xline(2062, HandleVisibility="off");
                xline(3093, HandleVisibility="off");
                xlabel("Times (ms)");
                ylabel("Power change (dB)");
                title("Power change from fixation baseline");
            
            end
            
            colorbar;
            sgtitle(obj.areainf + obj.condinflabel(tcond1) + "/" + obj.layeridlabel(layerid));
            legend;

        end

        function jLFPprobeINFO(obj, chx)

            imglfp1 = squeeze(mean(obj.x{1}(:, chx, :), 1));
        
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

        function flipObj = jVFLIP(obj, xch)

            data = permute(obj.x{1}(:, xch, :), [2, 3, 1]);
            flipObj = vFLIP2(data(:, :, :), 'DataType', 'raw_cut', 'fsample', 1000, 'intdist', 0.04, 'plot_result', true);
        
        end

    end

end


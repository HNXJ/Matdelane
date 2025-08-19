classdef jnwb < handle

    % JNWB: Passive omission task handler
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
        pgx2 = {};
        leakage = 0;
        freqlims = 0;
        overlap = 0;

        fmap = 0;
        tmap = 0;
        tbandlabels = ["Fix", "S1d1", "S2d2", "S3d3", "S4d4"];
        fbandlabels = ["Theta(2-7Hz)", "Alpha(8-12Hz)", "Beta(14-30Hz)", "GammaL(32-80Hz)", "GammaH(80+Hz)"];

        objectfile = "";
        channelinfo = cell(1);

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

            obj.tmapsignal = linspace(-obj.tpre, obj.tpost, size(obj.xm{1}, 3));

        end

        function jSave(obj, path, filename)

            filename = replace(filename, ["/", "\"], "");
            obj.objectfile = path + "\" + filename + ".mat";
            save(obj.objectfile, "obj", "-v7.3");

        end

        function jCleanLFPtrials(obj)

            for ik = 1:12

                dLFP = squeeze(mean(obj.x{ik}, 2));
                dCorrLFP = corr(dLFP');
                dCorrLFP = mean(dCorrLFP);
                dCorrLFP = (dCorrLFP - mean(dCorrLFP)) / std(dCorrLFP);
                dTRgood = abs(dCorrLFP) < 1.0;

                obj.x{ik} = obj.x{ik}(dTRgood, :, :);

            end

        end

        function jCleanMUAtrials(obj)

            for ik = 1:12

                dMUA = squeeze(mean(obj.xm{ik}, 2));
                dCorrMUA = corr(dMUA');
                dCorrMUA = mean(dCorrMUA);
                dCorrMUA = (dCorrMUA - mean(dCorrMUA)) / std(dCorrMUA);
                dTRgood = abs(dCorrMUA) < 1.0;

                obj.xm{ik} = obj.xm{ik}(dTRgood, :, :);

            end

        end

        function obj = jCalcTFRs(obj, channel_in_layer, tfrByChannel, tfrAllChannel) % TODO concat

            if ~exist("tfrByChannel", "var")

                tfrByChannel = 0;

            end

            if ~exist("tfrAllChannel", "var")

                tfrAllChannel = 0;

            end

            if tfrAllChannel

                [pgxt, xinfo] = jCalculateTFRallChannel(obj.x, channel_in_layer);

            else

                if tfrByChannel
    
                    [pgxt, xinfo] = jCalculateTFRbyChannel(obj.x, channel_in_layer);
    
                else
                
                    [pgxt, xinfo] = jCalculateTFR(obj.x, channel_in_layer);
    
                end

            end

            if sum(size(obj.pgx)) > 0

                obj.pgx2 = pgxt;

            else
            
                obj.pgx = pgxt;

            end

            obj.freqlims = xinfo.freqlims;
            obj.leakage = xinfo.leakage;

            obj.freqres = xinfo.freqres;
            obj.overlap = xinfo.overlap;
            obj.fmap = xinfo.fmap;

            obj.tmap = xinfo.tmap;
            obj.tbands = xinfo.tbands;
            obj.fbands = xinfo.fbands;

        end

        function obj = tmapsignalset(obj)

            obj.tmapsignal = linspace(-obj.tpre, obj.tpost, size(obj.x{1}, 3));

        end

        function jMUAplot(obj, condid, timewindow, smoothw, saveflag, channelsoi)

            if ~exist("timewindow", "var")

                obj.tmapsignalset();
                timewindow = [obj.tmapsignal(1), obj.tmapsignal(end)];

            end

            if ~exist("smoothw", "var")

                smoothw = 2;

            end

            if ~exist("saveflag", "var")

                saveflag = 0;

            end

            if ~exist("channelsoi", "var")

                channelsoi = 1:128;

            end

            figure;

            if length(condid) > 1

                clrx = obj.jMultiColors(length(condid));

                for condix = 1:length(condid)

                    imxm = squeeze(mean(obj.xm{condid(condix)}(:, channelsoi, :), 1));
                    imxmse = 2*squeeze(std(imxm)) / sqrt(size(imxm, 2));
                    imxm = squeeze(mean(imxm, 1));
                    imxm = (imxm - mean(imxm(250:500))) / std(imxm(500:2000));
                    imxm = smooth(imxm, smoothw);
                    imxmse = smooth(imxmse, smoothw);
                    plot(linspace(-500, 4250, 4750), imxm, "color", clrx{condix}, "DisplayName", obj.areainf + obj.condinflabel(condid(condix)));
                    
                    hold("on");
        
                    ximxmt = linspace(-500, 4250, 4750);
                    patch([ximxmt ximxmt(end:-1:1)], [imxmse + imxm; imxm(end:-1:1) - imxmse(end:-1:1)], clrx{condix}, "FaceAlpha",  .2, "HandleVisibility", "off");
        
                    xline(0, HandleVisibility="off");
                    xline(1031, HandleVisibility="off");
                    xline(2062, HandleVisibility="off");
                    xline(3093, HandleVisibility="off");
                    
                    title("MUAenv/Zsc/");
                    xlabel("Time (ms)");
                    ylabel("Z-score");
                    xlim(timewindow);

                end

                legend;

            else

                imxm = squeeze(mean(obj.xm{condid}(:, channelsoi, :), 1));
                imxmse = squeeze(std(imxm)) / sqrt(size(imxm, 1));
                imxm = squeeze(mean(imxm, 1));
                imxm = (imxm - mean(imxm(250:500))) / std(imxm(500:2000));
                imxm = smooth(imxm, smoothw);
                imxmse = smooth(imxmse, smoothw);
                plot(linspace(-500, 4250, 4750), imxm, "DisplayName", obj.areainf);
                
                hold("on");
    
                ximxmt = linspace(-500, 4250, 4750);
                patch([ximxmt ximxmt(end:-1:1)], [imxmse + imxm; imxm(end:-1:1) - imxmse(end:-1:1)], [.4 .5 .7], "FaceAlpha",  .2, "HandleVisibility", "off");
    
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

            if saveflag

                fname = replace(obj.areainf, "/", "") + "_mua_" + string(obj.condinflabel(condid(1)));
                print(gcf,'-vector','-dsvg', fname +".svg");

            end

        end

        function jMUAplotAll(obj, timewindow)

            if ~exist("timewindow", "var")

                obj.tmapsignalset();
                timewindow = [obj.tmapsignal(1), obj.tmapsignal(end)];

            end

            figure;
            
            for condid = 1:12

                imxm = squeeze(mean(obj.xm{condid}, 1));
                imxm = squeeze(mean(imxm, 1));
                imxm = (imxm - mean(imxm)) / std(imxm);
                subplot(3, 4, condid);
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

            sgt_ = char(obj.nwbFile);
            sgtitle(sgt_(6:end-4));

        end

        function jSUAplot(obj, condid, timewindow, units, trials)

            if ~exist("timewindow", "var")

                obj.tmapsignalset();
                timewindow = [obj.tmapsignal(1), obj.tmapsignal(end)];

            end

            if ~exist("units", "var")

                units = 1:size(obj.xs{condid}, 2);

            end

            if ~exist("trials", "var")

                trials = 1:size(obj.xs{condid}, 1);

            end

            figure;
            
            imxm = squeeze(mean(obj.xs{condid}(trials, units, :), 1));
            imxm = squeeze(mean(imxm, 1))*1000;
            imxm = smooth(imxm, 100, "sgolay");
            % imxm = imxm + randn(size(imxm))*10;
            % imxm = (imxm - mean(imxm)) / std(imxm);
            plot(linspace(-500, 4250, 4750), imxm, "DisplayName", obj.areainf);
            
            hold("on");
            xline(0, HandleVisibility="off");
            xline(1031, HandleVisibility="off");
            xline(2062, HandleVisibility="off");
            xline(3093, HandleVisibility="off");
            
            title("SUAenv/est" + obj.condinflabel(condid));
            xlabel("Time (ms)");
            ylabel("Sp/s");
            xlim(timewindow);
            
            legend;

        end

        function jRastrogram(obj, condid, timewindow, units, trials)

            if ~exist("timewindow", "var")

                obj.tmapsignalset();
                timewindow = [obj.tmapsignal(1), obj.tmapsignal(end)];

            end

            if ~exist("units", "var")

                units = 1:size(obj.xs{condid}, 2);

            end

            if ~exist("trials", "var")

                trials = 1:size(obj.xs{condid}, 1);

            end

            figure;
            
            imxm = squeeze(mean(obj.xs{condid}(trials, units, :), 1));
            imagesc(imxm, "XData", linspace(-500, 4250, 4750));
            colormap("gray");
            
            hold("on");
            xline(0, HandleVisibility="off", LineWidth=1.0, Color=[1 1 1]);
            xline(1031, HandleVisibility="off", LineWidth=1.0, Color=[1 1 1]);
            xline(2062, HandleVisibility="off", LineWidth=1.0, Color=[1 1 1]);
            xline(3093, HandleVisibility="off", LineWidth=1.0, Color=[1 1 1]);
            
            title("Rastrogram/" + obj.condinflabel(condid));
            xlabel("Time (ms)");
            ylabel("Neuron id");
            xlim(timewindow);

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
            % tfr1e = squeeze(std(obj.pgx{tcond1, layerid}));
            
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
            cls = zeros(5, 3);
            cls(1, :) = [0 0 1];
            cls(2, :) = [1 0 0];
            cls(3, :) = [1 0.5 0];
            cls(4, :) = [1 0 1];
            cls(5, :) = [0 0.8 0.4];
            
            for fband = 1:5
            
                tfr1 = squeeze(mean(obj.pgx{tcond1, layerid}(:, obj.fbands{fband}, :), 1));
                tfr1 = squeeze(mean(tfr1, 1));
                tfr1 = smooth(tfr1, 10);
                baselinecrx = mean(tfr1(tbaseline));

                tfr1 = tfr1 / baselinecrx;

                ne1 = length(obj.fbands{fband}); % spectral res. bands
                ne2 = size(obj.pgx{tcond1, layerid}, 1); % trial count
                ne3 = size(obj.x{tcond1}, 2); % channel count

                tfr1e = squeeze(std(obj.pgx{tcond1, layerid}(:, obj.fbands{fband}, :) / baselinecrx));
                tfr1e = squeeze(mean(tfr1e, 1)) / sqrt(ne1 + ne2 + ne3);
                tfr1e = smooth(tfr1e, 10);

                cl = cls(fband, :);
                y1s = 10*log10(tfr1);
                plot(obj.tmap, y1s, "DisplayName", obj.fbandlabels(fband), "LineWidth", 0.5, "Color", cl);
                hold("on");

                stx = tfr1 + tfr1e;
                sty = tfr1 - tfr1e;
                stx(stx <= 0) = 1e-2;
                sty(sty <= 0) = 1e-2;

                stx = 10*log10(stx);
                sty = 10*log10(sty);

                patch([obj.tmap; obj.tmap(end:-1:1)], [stx; sty(end:-1:1)], cl, "EdgeColor", "none", "FaceColor", cl, "FaceAlpha", 0.5, "HandleVisibility", "off");

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
            sgtitle(obj.areainf + obj.condinflabel(tcond1) + "/" + obj.layeridlabel(layerid));
            legend;

        end

        function jLFPprobeINFO(obj, chx, condx, txlims)

            if ~exist("condx", "var")

                condx = 1;

            end

            if ~exist("txlims", "var")

                txlims = [0, obj.tpost];

            end

            imglfp1 = squeeze(mean(obj.x{condx}(:, chx, :), 1));
        
            figure;
            subplot(4, 1, 1);
            imagesc(imglfp1, "YData", chx);
            xlim(txlims);
            % clim([-2.5 2.5]);
            xlabel("Time (ms)");
            ylabel("Channel");
            title("LFP average by trial");
        
            imglfp2 = smoothdata2(imglfp1, "gaussian", {10, 20});
            imgcsd1 = (-imglfp2(1:end-2, :) + 2*imglfp2(2:end-1, :) - imglfp2(3:end, :));
            imgcsd1 =  smoothdata2(imgcsd1, "gaussian", {7, 20});

            % imgcsd1 = (imgcsd1);
        
            subplot(4, 1, 2);
            imagesc(imgcsd1, "YData", chx);
            xlim(txlims);
            clim([-2.5 2.5]);
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

        function jLFPmuaINFO(obj, chx, condx, txlims)

            if ~exist("condx", "var")

                condx = 1;

            end

            if ~exist("txlims", "var")

                txlims = [0, obj.tpost];

            end

            imglfp1 = squeeze(mean(obj.x{condx}(:, chx, :), 1));
            imgmua1 = squeeze(mean(obj.xm{condx}(:, chx, :), 1));
        
            figure;
            subplot(4, 1, 1);
            imagesc(imglfp1, "YData", chx);
            xlim(txlims);
            % clim([-2.5 2.5]);
            xlabel("Time (ms)");
            ylabel("Channel");
            title("LFP average by trial");

            subplot(4, 1, 2);
            imagesc(imgmua1, "YData", chx);
            xlim(txlims);
            clim([-2.5 2.5]);
            xlabel("Time (ms)");
            ylabel("Channel");%clim([-1 1]);
            title("MUA average by trial");
        
            corrx = corr(imglfp1', imgmua1', "Type", "Spearman");
            subplot(2, 2, 3);
            imagesc(corrx);
            xlabel("Channel");
            ylabel("Channel");
            title("LFP-MUA corr average by trial");
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

        function [expvars, layerinf] = jCalcPEV(obj, layerid, condinf)
            
            layerinf = obj.layeridlabel(layerid) + " layer";
            expvars = cell(1, 5);
            
            for fband = 1:5
            
                x1 = obj.pgx{condinf(1), layerid}(:, obj.fbands{fband}, :);
                x2 = obj.pgx{condinf(2), layerid}(:, obj.fbands{fband}, :);
                
                [N1, nF, nT] = size(x1);
                N2 = size(x2, 1);
                
                data = zeros(N1+N2, nF, nT);
                data(1:N1, :, :) = x1;
                data(N1+1:N1+N2, :, :) = x2;
                % data = jSmooth(data, 50);
                
                groupIDs = [ones(1, N1), ones(1, N2)*2];
                [expv, ~, ~, ~, ~] = jPEV(data, groupIDs, 1);
                expvars{fband} = squeeze(expv)*100;
            
            end

        end

        function [expvars, layerinf] = jCalcPEVs(obj, layerid, condinf)
            
            layerinf = obj.layeridlabel(layerid) + " layer";
            expvars = cell(1, 5);
            
            for fband = 1:5
            
                x1 = obj.pgx{condinf(1), layerid}(:, obj.fbands{fband}, :);
                x2 = obj.pgx{condinf(2), layerid}(:, obj.fbands{fband}, :);
                x3 = obj.pgx{condinf(3), layerid}(:, obj.fbands{fband}, :);
                
                [N1, nF, nT] = size(x1);
                N2 = size(x2, 1);
                N3 = size(x3, 1);
                
                data = zeros(N1+N2, nF, nT);
                data(1:N1, :, :) = x1;
                data(N1+1:N1+N2, :, :) = x2;
                data(N1+N2+1:N1+N2+N3, :, :) = x3;
                % data = jSmooth(data, 50);
                
                groupIDs = [ones(1, N1), ones(1, N2)*2, ones(1, N3)*3];
                [expv, ~, ~, ~, ~] = jPEV(data, groupIDs, 1);
                expvars{fband} = squeeze(expv)*100;
            
            end

        end

        function jPEVplot(obj, expvars, layerinf, condinf)

            figure;
            
            for fband = 1:5
            
                subplot(5, 1, fband);
                yt = mean(expvars{fband}, 1);
                yt = yt - min(yt);
                st = std(expvars{fband}) / sqrt(size(expvars{fband}, 1));
                plot(obj.tmap, yt);hold("on");
                stx = smooth(yt + 2*st, 2);stx(stx<0) = 0;
                sty = smooth(yt - 2*st, 2);sty(sty<0) = 0;
                cl = [0.9 0.7 0.7];
                plot(obj.tmap, stx, "Color", cl);
                plot(obj.tmap, sty, "Color", cl);
                patch([obj.tmap', obj.tmap(end:-1:1)'], [stx;sty(end:-1:1)], cl);
                xline(0);
                xline(1031);
                xline(2062);
                xline(3093);
                xlabel("Time(ms)");ylabel("PEV(%)");
                title(obj.fbandlabels(fband));
            
            end
            
            sgtitle("Area:" + obj.areainf + " " + obj.condinflabel(condinf(1)) + "-vs-" + obj.condinflabel(condinf(2)) + " PEV/TFR/+-2SEM/fRes=" + num2str(obj.freqres) + "Hz/ovlrp=." + num2str(obj.overlap) + " " + layerinf);

        end

        function flipObj = jVFLIP(obj, xch, txlim, cond, dist, plotflag)

            if ~exist('cond', 'var')

                cond = 1;

            end

            if ~exist('dist', 'var')

                dist = 0.040;

            end

            if ~exist('plotflag', 'var')

                plotflag = 1;

            end

            if ~exist('txlim', 'var')

                data = permute(obj.x{cond}(:, xch, :), [2, 3, 1]);

            else

                data = permute(obj.x{cond}(:, xch, txlim), [2, 3, 1]);

            end
            
            flipObj = vFLIP2(data(:, :, :), 'DataType', 'raw_cut', 'fsample', 1000, 'intdist', dist, 'plot_result', plotflag);
        
        end

        function [tfr1x, tfr2x] = jTFR2dScatters(obj, tcond1, tcond2, layerid, tbaseline, txlims) % TODO

            if ~exist("tbaseline", "var")
            
                tbaseline = obj.tbands{1};

            end

            if ~exist("txlims", "var")
            
                txlims = [obj.tmap(1) obj.tmap(end)];

            end
           
            tfr1 = squeeze(mean(obj.pgx{tcond1, layerid}, 1));
            tfr2 = squeeze(mean(obj.pgx{tcond2, layerid}, 1));
            
            for ik = 1:size(tfr1, 1)
            
                tfr1(ik, :) = tfr1(ik, :) / mean(tfr1(ik, tbaseline));
            
            end
                        
            for ik = 1:size(tfr2, 1)
            
                tfr2(ik, :) = tfr2(ik, :) / mean(tfr2(ik, tbaseline));
            
            end

            tfr1x = zeros([5, size(tfr1, 2)]);
            tfr2x = zeros([5, size(tfr2, 2)]);

            cls = zeros(5, 3);
            cls(1, :) = [0 0 1];
            cls(2, :) = [1 0 0];
            cls(3, :) = [1 0.5 0];
            cls(4, :) = [1 0 1];
            cls(5, :) = [0 0.8 0.4];
            
            for fband = 1:5
            
                tfr1t = squeeze(mean(obj.pgx{tcond1, layerid}(:, obj.fbands{fband}, :), 1));
                tfr1t = squeeze(mean(tfr1t, 1));
                tfr1t = smooth(tfr1t, 10);
                baselinecrx = mean(tfr1t(tbaseline));

                tfr1t = tfr1t / baselinecrx;
                tfr1x(fband, :) = tfr1t;

                tfr2t = squeeze(mean(obj.pgx{tcond2, layerid}(:, obj.fbands{fband}, :), 1));
                tfr2t = squeeze(mean(tfr2t, 1));
                tfr2t = smooth(tfr2t, 10);
                baselinecrx = mean(tfr2t(tbaseline));

                tfr2t = tfr2t / baselinecrx;
                tfr2x(fband, :) = tfr2t;

            end
           
        end

        function clrxs = jMultiColors(~, N)

            clrxs = cell(1, N);

            for ik = 1:N

                clrxs{ik} = [mod(ik, 2) mod(floor(ik/2), 2), mod(floor(ik/4), 2)]/2 + 0.5;

            end

        end

    end

end


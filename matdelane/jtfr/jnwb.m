classdef jnwb

    % JNWB: Time-Frequency analysis on signals
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

        xs = {};
        cs = {};
        
    end
    
    methods

        function obj = jnwb(nwbFile, areainf, tpre, tpost, probeid, convFlag)

            nwb = nwbRead(nwbFile);
            [obj.c, obj.x] = jOGLOSignals(nwb, "omission_glo_passive", tpre, tpost, probeid, "lfp");
            [obj.cm, obj.xm] = jOGLOSignals(nwb, "omission_glo_passive", tpre, tpost, probeid, "muae");
            [obj.cs, obj.xs] = jOGLOUnits(nwb, "omission_glo_passive", tpre, tpost, convFlag);

            disp(obj.c{1}.session);
            obj.areainf = areainf;

        end
       
        function pgx = jTFR(obj, x)

            pgx = x;

        end

        function jTFRplot(obj, pgx, tcond1, layerid)

            tcond1 = 11;
            layerid = 4;
            tbaseline3 = obj.tbands{3}(end-15:end-10);
            txlims = [obj.tmap(obj.tbands{3}(1)) obj.tmap(obj.tbands{4}(end))];
            
            figure;
            subplot(2, 1, 1);
            tfr1 = squeeze(mean(pgx{tcond1, layerid}, 1));
            
            for ik = 1:size(tfr1, 1)
            
                tfr1(ik, :) = tfr1(ik, :) / mean(tfr1(ik, tbaseline3));
            
            end
            
            imagesc(10*log10(tfr1), "XData", tmap, "YData", fmap);
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
            
                yline(fmap(fbands{fband}(1)), "Color", [1 0 0]);
            
            end
            
            xlabel("Time (ms)");
            subplot(2, 1, 2);
            
            for fband = 1:5
            
                tfr1 = squeeze(mean(pgx{tcond1, layerid}(:, fbands{fband}, :), 1));
                tfr1 = squeeze(mean(tfr1, 1));
                tfr1 = smooth(tfr1, 20);
                tfr1 = tfr1 / mean(tfr1(tbaseline3));
                plot(tmap, 10*log10(tfr1), "DisplayName", fbandlabels(fband), "LineWidth", 2);
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
            sgtitle(areainf + condinflabel(tcond1) + "/" + layeridlabel(layerid));
            legend;

        end

    end
end


function [probeInfo, signalList] = jNWBSignals(nwb, probe, task, t_pre_ms, t_post_ms, signalMode)

    if ~exist('probe', 'var')

        probe = 0;

    end

    if ~exist('t_pre_ms', 'var')

        t_pre_ms = 1000;

    end

    if ~exist('t_post_ms', 'var')

        t_post_ms = 4000;

    end

    if ~exist('mode', 'var')

        signalMode = "lfp";

    end

    if ~exist('task', 'var')

        tasks = nwb.intervals.keys();
        task_warning = "->Enter an specific task name from this file. \n-->This file contains these tasks : ";
        for i = 1:numel(tasks)
            task_warning = task_warning + " \n--->" + string(tasks{i});
        end

        wtx = sprintf(task_warning);
        warning(wtx);
        signalList = {};
        probeInfo = "Null";
        return;

    else

        task = string(task);

        try
        
            correct = nwb.intervals.get(task);

        catch

            tasks = nwb.intervals.keys();
            task_warning = "->Enter an specific task name from this file. \n-->This file contains these tasks : ";
            for i = 1:numel(tasks)
                task_warning = task_warning + " \n--->" + string(tasks{i});
            end
    
            wtx = sprintf(task_warning);
            warning(wtx);
            signalList = {};
            probeInfo = "Null";
            return;

        end

    end

    signalList = cell(1, 1);
    jNWBTaskDetails(nwb, task);

    probeLabel = ['probe', char(65 + probe)];
    areaName = nwb.general_extracellular_ephys.get(probeLabel).location{1};
    disp(['Area(s): ', areaName]);
    probeInfo = areaName;

    blocks = nwb.intervals.get(task).vectordata.get("task_block_number").data(:);
    correct = nwb.intervals.get(task).vectordata.get("correct").data(:);
    conditions = nwb.intervals.get(task).vectordata.get("task_condition_number").data(:);
    stims = nwb.intervals.get(task).vectordata.get("stimulus_number").data(:);

    conditionlist = unique(conditions);

    try
    
        sig = nwb.acquisition.get("probe_" + num2str(probe) + "_" + signalMode).electricalseries.get("probe_" + num2str(probe) + "_" + signalMode + "_data").data;

    catch

        disp("->MATNWB: Problem with this file while loading the signal mode : " + signalMode);
        signalList = {};
        return;

    end
    
    Nchannel = sig.dims(1);
    stime = nwb.intervals.get(task).start_time.data(:);
    stimeind = floor(stime*1000);

    for condition = conditionlist'
            
        b = (conditions == condition) & (correct == 1) & (stims == 1);
        b = find(b);
        temp_signals = zeros(numel(b), Nchannel, t_pre_ms + t_post_ms);
        cnt = 0;

        for i = b'

            cnt = cnt + 1;
            temp_signals(cnt, :, :) = sig(:, stimeind(i) - t_pre_ms + 1:stimeind(i) + t_post_ms);

        end

        signalList{condition} = temp_signals;

    end

    disp("_");

end
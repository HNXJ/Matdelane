function [signalListL, signalList] = jOGLOSignals(nwb, task, t_pre_ms, t_post_ms, probeId, signalMode)

    if ~exist('t_pre_ms', 'var')

        t_pre_ms = 1000;

    end

    if ~exist('t_post_ms', 'var')

        t_post_ms = 4000;

    end

    if ~exist('probeId', 'var')

        probeId = 0;

    end

    if ~exist('signalMode', 'var')

        signalMode = "lfp";

    end

    disp("->Signal mode : " + signalMode);

    if ~exist('task', 'var')

        tasks = nwb.intervals.keys();
        task_warning = "->Enter an specific task name from this file. \n-->This file contains these tasks : ";
        for i = 1:numel(tasks)
            task_warning = task_warning + " \n--->" + string(tasks{i});
        end

        warning(task_warning);
        signalList = {};
        return;

    else

        task = string(task);

    end

    signalListL = cell(1, 1);

    try

        rfdata = jReceptiveFieldIdx(nwb, probeId + 1);

    catch

        warning("->No RF-task exists in this NWB file.");
        rfdata = zeros([128, 4]);

    end

    areaN = 0;

    probeLabel = ['probe', char(65 + probeId)];
    areaName = nwb.general_extracellular_ephys.get(probeLabel).location{1};
    areas = strsplit(areaName, ',');
    areaCount = length(areas);

    channelW = 128/areaCount;
    
    for j = 1:areaCount

        k = areas{j};
        area = k(~isspace(k));
        areaN = areaN + 1;  

        signalListL{areaN}.ids = ceil(((j-1)*channelW)+1:j*channelW);
        signalListL{areaN}.name = area;      
        signalListL{areaN}.session = nwb.identifier;
        signalListL{areaN}.rfdata = rfdata(signalListL{areaN}.ids, :);

    end

    correct = nwb.intervals.get(task).vectordata.get("correct").data(:);
    conditions = nwb.intervals.get(task).vectordata.get("task_condition_number").data(:);
    stims = nwb.intervals.get(task).vectordata.get("stimulus_number").data(:);

    try
    
        sig = nwb.acquisition.get("probe_" + num2str(probeId) + "_" + signalMode).electricalseries.get("probe_" + num2str(probeId) + "_" + signalMode + "_data").data;

    catch

        disp("->LFP not found.");
        signalList = {};
        return;

    end
    
    stime = nwb.intervals.get(task).start_time.data(:);
    stimeind = floor(stime*1000);
    disp(num2str(areaN) + " areas identified.");

    conditionlist = [0, 2, 3, 4, 5, 7, 8, 9, 10, 26, 34, 42, 50];
    conditionlist_rrxr = 35:2:42;
    conditionlist_rrrx = [36:2:42, 43:50];
    signalList = cell(1, numel(conditionlist) - 1);

    parfor k = 2:numel(conditionlist)
        
        if k < 12
            b = (conditions <= conditionlist(k)) & (conditions > conditionlist(k-1)) & (correct == 1) & (stims == 2);
        elseif k == 12
            b = (ismember(conditions, conditionlist_rrxr)) & (correct == 1) & (stims == 2);
        else
            b = (ismember(conditions, conditionlist_rrrx)) & (correct == 1) & (stims == 2);
        end

        b = find(b);
        temp_signals = zeros(numel(b), 128, t_pre_ms + t_post_ms);
        cnt = 0;

        for i = b'

            cnt = cnt + 1;
            temp_signals(cnt, :, :) = sig(1:128, stimeind(i) - t_pre_ms + 1:stimeind(i) + t_post_ms);

        end

        signalList{k-1} = temp_signals;
        fprintf("%d ", k-1);

    end

    fprintf("\nDone.\n");

end
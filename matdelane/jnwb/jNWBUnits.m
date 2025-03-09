function [signalListL, signalList] = jNWBUnits(nwb, task, t_pre_ms, t_post_ms)

    if ~exist('t_pre_ms', 'var')

        t_pre_ms = 1000;

    end

    if ~exist('t_post_ms', 'var')

        t_post_ms = 4000;

    end

    if ~exist('task', 'var')

        tasks = nwb.intervals.keys();
        task_warning = "->Enter an specific task name from this file. \n-->This file contains these tasks : ";
        for i = 1:numel(tasks)
            task_warning = task_warning + " \n--->" + string(tasks{i});
        end

        warning(task_warning);
        signalListL = {};
        signalList = {};
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
            signalListL = {};
            return;

        end

    end

    signalListL = cell(1, 1);
    probeCount = length(nwb.general_extracellular_ephys.keys());
    
    areaN = 0;
    goodUnits = find(nwb.units.vectordata.get("quality").data(:) == 1);
    chidUnits = nwb.units.vectordata.get("peak_channel_id").data(goodUnits);

    for i = 1:probeCount
        
        probeLabel = ['probe', char(65 + i-1)];
        areaName = nwb.general_extracellular_ephys.get(probeLabel).location{1};
        areas = strsplit(areaName, ',');
        areaCount = length(areas);

        channelW = 128/areaCount;
        
        for j = 1:areaCount

            k = areas{j};
            area = k(~isspace(k));
            channelL = ceil(channelW*(j-1)) + (i-1)*128 + 1;
            channelR = ceil(channelW*j) + (i-1)*128;
           
            areaN = areaN + 1;  
            signalListL{areaN}.ids = find(chidUnits > channelL & chidUnits < channelR);
            signalListL{areaN}.name = area;      

        end

    end

    correct = nwb.intervals.get(task).vectordata.get("correct").data(:);
    conditions = nwb.intervals.get(task).vectordata.get("task_condition_number").data(:);
    stims = nwb.intervals.get(task).vectordata.get("stimulus_number").data(:);
    conditionlist = unique(conditions);

    try
    
        sig = nwb.processing.get("convolved_spike_train").nwbdatainterface.get("convolved_spike_train_data").data;

    catch

        disp("->Kilosort2.0: Single unit neurons were not detected in this probe.");
        return;

    end
    
    stime = nwb.intervals.get(task).start_time.data(:);
    stimeind = floor(stime*1000);
    disp(num2str(areaN) + " areas identified.");
    signalList = cell(1, numel(conditionlist));

    NgoodUnitsID = goodUnits;
    NgoodUnits = numel(NgoodUnitsID);

    parfor condition = conditionlist'
        
        b = (conditions == condition) & (correct == 1) & (stims == 2);
        b = find(b);
        temp_signals = zeros(numel(b), NgoodUnits, t_pre_ms + t_post_ms);
        cnt = 0;

        for i = b'

            cnt = cnt + 1;
            temp_signals(cnt, :, :) = sig(NgoodUnitsID, stimeind(i) - t_pre_ms + 1:stimeind(i) + t_post_ms);

        end

        signalList{condition} = temp_signals;
        fprintf("%d ", condition);

    end

    fprintf("\nDone.\n");

end
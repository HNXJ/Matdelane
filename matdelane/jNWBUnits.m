function signalList = jNWBUnits(nwb, probe, task, t_pre_ms, t_post_ms)

    if ~exist('probe', 'var')

        probe = 0;

    end

    if ~exist('t_pre_ms', 'var')

        t_pre_ms = 1000;

    end

    if ~exist('t_post_ms', 'var')

        t_pre_ms = 4000;

    end

    if ~exist('task', 'var')

        tasks = nwb.intervals.keys();
        task_warning = "->Enter an specific task name from this file. \n-->This file contains these tasks : ";
        for i = 1:numel(tasks)
            task_warning = task_warning + " \n--->" + string(tasks{i});
        end

        warning(sprintf(task_warning));
        signalList = {};
        return;

    else

        task = string(task);

    end

    signalList = cell(1, 1);
    jNWBTaskDetails(nwb, task);

    probeLabel = ['probe', char(65 + probe)];
    areaName = nwb.general_extracellular_ephys.get(probeLabel).location{1};
    disp(['Area(s): ', areaName]);

    blocks = nwb.intervals.get(task).vectordata.get("task_block_number").data(:);
    correct = nwb.intervals.get(task).vectordata.get("correct").data(:);
    conditions = nwb.intervals.get(task).vectordata.get("task_condition_number").data(:);
    stims = nwb.intervals.get(task).vectordata.get("stimulus_number").data(:);

    blocklist = unique(blocks);
    conditionlist = unique(conditions);
    
    N = size(blocks, 1);

    try
    
        sig = nwb.processing.get("convolved_spike_train").nwbdatainterface.get("convolved_spike_train_data").data;

    catch

        disp("->Kilosort2.0: Single unit neurons were not detected in this probe.");
        signalList = {};
        return;

    end
    
    goodUnits = nwb.units.vectordata.get("quality").data(:);
    chidUnits = nwb.units.vectordata.get("peak_channel_id").data(:) >= probe*128 & nwb.units.vectordata.get("peak_channel_id").data(:) < (probe+1)*128;
    stime = nwb.intervals.get(task).start_time.data(:);
    stimeind = floor(stime*1000);
        
    goodUnits = (goodUnits & chidUnits);
    NgoodUnits = sum(goodUnits);
    NgoodUnitsID = find(goodUnits == 1);

    for block = blocklist'
    
        for condition = conditionlist'
            
            b = (blocks == block) & (conditions == condition) & (correct == 1) & (stims == 1);
            b = find(b);
            temp_signals = zeros(numel(b), NgoodUnits, t_pre_ms + t_post_ms);
            cnt = 0;

            for i = b'

                cnt = cnt + 1;
                temp_signals(cnt, :, :) = sig(NgoodUnitsID, stimeind(i) - t_pre_ms + 1:stimeind(i) + t_post_ms);

            end

            signalList{block, condition} = temp_signals;

        end
    
    end

    disp("_");

end
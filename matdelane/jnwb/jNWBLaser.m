function signalList = jNWBLaser(nwb, laser, task, t_pre_ms, t_post_ms)

    if ~exist('laser', 'var')

        laser = 0;

    end

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

        warning(sprintf(task_warning));
        signalList = {};
        return;

    else

        task = string(task);

    end

    signalList = cell(1, 1);
    jNWBTaskDetails(nwb, task);

    blocks = nwb.intervals.get(task).vectordata.get("task_block_number").data(:);
    correct = nwb.intervals.get(task).vectordata.get("correct").data(:);
    conditions = nwb.intervals.get(task).vectordata.get("task_condition_number").data(:);
    stims = nwb.intervals.get(task).vectordata.get("stimulus_number").data(:);

    blocklist = unique(blocks);
    conditionlist = unique(conditions);
    
    N = size(blocks, 1);

    try
    
        sig = nwb.acquisition.get("laser_" + num2str(laser) + "_tracking").timeseries.get("laser_" + num2str(laser) + "_tracking_data").data;

    catch

        disp("->OPTRODE: Specified laser " + num2str(laser) + " not found in this file.");
        signalList = {};
        return;

    end
    
    stime = nwb.intervals.get(task).start_time.data(:);
    stimeind = floor(stime*1000);

    for block = blocklist'
    
        for condition = conditionlist'
            
            b = (blocks == block) & (conditions == condition) & (correct == 1) & (stims == 2);
            b = find(b);
            temp_signals = zeros(numel(b), t_pre_ms + t_post_ms);
            cnt = 0;

            for i = b'

                cnt = cnt + 1;
                temp_signals(cnt, :) = sig(stimeind(i) - t_pre_ms + 1:stimeind(i) + t_post_ms);

            end

            signalList{block, condition} = temp_signals;

        end
    
    end

    disp("_");

end
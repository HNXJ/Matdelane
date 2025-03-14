function jNWBTaskDetails(nwb, taskname)

    if ~exist('taskname', 'var')
    
        tasks = nwb.intervals.keys();
        task_warning = "->Enter an specific task name from this file. \n-->This file contains these tasks : ";
        for i = 1:numel(tasks)
            task_warning = task_warning + " \n--->" + string(tasks{i});
        end

        wtx = sprintf(task_warning);
        warning(wtx);
        signalList = {};
        return;

    else

        taskname = string(taskname);
        
        try
        
            correct = nwb.intervals.get(taskname);

        catch

            tasks = nwb.intervals.keys();
            task_warning = "->Enter an specific task name from this file. \n-->This file contains these tasks : ";
            for i = 1:numel(tasks)
                task_warning = task_warning + " \n--->" + string(tasks{i});
            end
    
            wtx = sprintf(task_warning);
            warning(wtx);
            signalList = {};
            return;

        end
    
    end

    info_text = "";
    
    correct = nwb.intervals.get(taskname).vectordata.get("correct").data(:);
    blocks = nwb.intervals.get(taskname).vectordata.get("task_block_number").data(:);
    conditions = nwb.intervals.get(taskname).vectordata.get("task_condition_number").data(:);
    stims = nwb.intervals.get(taskname).vectordata.get("stimulus_number").data(:);

    x1 = (correct == 1) & (stims == 1);
    x2 = numel(unique(blocks));
    x3 = numel(unique(conditions));
    x4 = (stims == 1);

    info_text = info_text + sprintf("\n>" + taskname + ":\n->Total correct trials: %d\n", sum(x1));
    info_text = info_text + sprintf("\n-->Total blocks: %d\n", sum(x2));
    info_text = info_text + sprintf("\n-->Total conditions: %d\n", sum(x3));

    for i = 1:max(conditions)

        cond_cnt = sum((conditions == i) & (stims == 1));
        cond_correct = sum((conditions == i) & (correct == 1) & (stims == 1));
        info_text = info_text + sprintf("\n--->Condition %d : %d trials, %d correct \n", i, sum(cond_cnt), sum(cond_correct));

    end

    info_text = info_text + sprintf("\n->Total trials: %d\n", sum(x4));
    disp(info_text);

end
function signalList = jOGLOUnits(nwb, probe)

    if ~exist('probe', 'var')

        probe = 0;

    end

    % jOGLOTaskDetails(nwb);

    probeLabel = ['probe', char(65 + probe)];

    try

        areaName = nwb.general_extracellular_ephys.get(probeLabel).location{1};

    catch

        signalList = struct();
        signalList.area = 'null';
        return;

    end

    disp(['Area(s): ', areaName]);

    blocks = nwb.intervals.get("omission_glo_passive").vectordata.get("task_block_number").data(:);
    correct = nwb.intervals.get("omission_glo_passive").vectordata.get("correct").data(:);
    omission = nwb.intervals.get("omission_glo_passive").vectordata.get("is_omission").data(:);
    stims = nwb.intervals.get("omission_glo_passive").vectordata.get("stimulus_number").data(:);
    
    N = size(blocks, 1);

    try
    
        sig = nwb.processing.get("convolved_spike_train").nwbdatainterface.get("convolved_spike_train_data").data;

    catch

        disp("->Kilosort2.0: Single unit neurons were not detected in this probe.");
        signalList = -1;
        return;

    end
    
    goodUnits = nwb.units.vectordata.get("quality").data(:);
    chidUnits = nwb.units.vectordata.get("peak_channel_id").data(:) >= probe*128 & nwb.units.vectordata.get("peak_channel_id").data(:) < (probe+1)*128;
    stime = nwb.intervals.get("omission_glo_passive").start_time.data(:);
    stimeind = floor(stime*1000);
    
    b2l = (blocks == 2) & (correct == 1) & (isnan(omission)) & (stims == 3);
    b2o2 = (blocks == 2) & (correct == 1) & (omission == 1) & (stims == 3);
    b2o3 = (blocks == 2) & (correct == 1) & (omission == 1) & (stims == 4);
    b2o4 = (blocks == 2) & (correct == 1) & (omission == 1) & (stims == 5);
    
    b4l = (blocks == 4) & (correct == 1) & (isnan(omission)) & (stims == 3);
    b4o2 = (blocks == 4) & (correct == 1) & (omission == 1) & (stims == 3);
    b4o3 = (blocks == 4) & (correct == 1) & (omission == 1) & (stims == 4);
    b4o4 = (blocks == 4) & (correct == 1) & (omission == 1) & (stims == 5);
    
    b5l = (blocks == 5) & (correct == 1) & (isnan(omission)) & (stims == 3);
    b5o2 = (blocks == 5) & (correct == 1) & (omission == 1) & (stims == 3);
    b5o3 = (blocks == 5) & (correct == 1) & (omission == 1) & (stims == 4);
    b5o4 = (blocks == 5) & (correct == 1) & (omission == 1) & (stims == 5);
    
    goodUnits = (goodUnits & chidUnits);
    NgoodUnits = sum(goodUnits);
    NgoodUnitsID = find(goodUnits == 1);

    sigsb2l = zeros(sum(b2l), NgoodUnits, 6000);
    sigsb2o2 = zeros(sum(b2o2), NgoodUnits, 6000);
    sigsb2o3 = zeros(sum(b2o3), NgoodUnits, 6000);
    sigsb2o4 = zeros(sum(b2o4), NgoodUnits, 6000);
    
    sigsb4l = zeros(sum(b4l), NgoodUnits, 6000);
    sigsb4o2 = zeros(sum(b4o2), NgoodUnits, 6000);
    sigsb4o3 = zeros(sum(b4o3), NgoodUnits, 6000);
    sigsb4o4 = zeros(sum(b4o4), NgoodUnits, 6000);
    
    sigsb5l = zeros(sum(b5l), NgoodUnits, 6000);
    sigsb5o2 = zeros(sum(b5o2), NgoodUnits, 6000);
    sigsb5o3 = zeros(sum(b5o3), NgoodUnits, 6000);
    sigsb5o4 = zeros(sum(b5o4), NgoodUnits, 6000);
    
    cntl = 0;
    cnto2 = 0;
    cnto3 = 0;
    cnto4 = 0;
    
    for i = 1:N
    
        if b2l(i)
    
            cntl = cntl + 1;
            sigsb2l(cntl, :, :) = sig(NgoodUnitsID, stimeind(i-1) - 1000 + 1:stimeind(i-1) + 5000);
    
        elseif b2o2(i)
    
            cnto2 = cnto2 + 1;
            sigsb2o2(cnto2, :, :) = sig(NgoodUnitsID, stimeind(i-1) - 1000 + 1:stimeind(i-1) + 5000);
    
        elseif b2o3(i)
    
            cnto3 = cnto3 + 1;
            sigsb2o3(cnto3, :, :) = sig(NgoodUnitsID, stimeind(i-2) - 1000 + 1:stimeind(i-2) + 5000);
    
        elseif b2o4(i)
    
            cnto4 = cnto4 + 1;
            sigsb2o4(cnto4, :, :) = sig(NgoodUnitsID, stimeind(i-3) - 1000 + 1:stimeind(i-3) + 5000);
    
        end
    
        if mod(i, 2000) == 1
            fprintf(num2str(i) + "-");
        end
    
    end
    
    disp("_");

    cntl = 0;
    cnto2 = 0;
    cnto3 = 0;
    cnto4 = 0;
    
    for i = 1:N
    
        if b4l(i)
    
            cntl = cntl + 1;
            sigsb4l(cntl, :, :) = sig(NgoodUnitsID, stimeind(i-1) - 1000 + 1:stimeind(i-1) + 5000);
    
        elseif b4o2(i)
    
            cnto2 = cnto2 + 1;
            sigsb4o2(cnto2, :, :) = sig(NgoodUnitsID, stimeind(i-1) - 1000 + 1:stimeind(i-1) + 5000);
    
        elseif b4o3(i)
    
            cnto3 = cnto3 + 1;
            sigsb4o3(cnto3, :, :) = sig(NgoodUnitsID, stimeind(i-2) - 1000 + 1:stimeind(i-2) + 5000);
    
        elseif b4o4(i)
    
            cnto4 = cnto4 + 1;
            sigsb4o4(cnto4, :, :) = sig(NgoodUnitsID, stimeind(i-3) - 1000 + 1:stimeind(i-3) + 5000);
    
        end
    
        if mod(i, 2000) == 1
            fprintf(num2str(i) + "-");
        end
    
    end

    disp("_");

    cntl = 0;
    cnto2 = 0;
    cnto3 = 0;
    cnto4 = 0;
    
    for i = 1:N
    
        if b5l(i)
    
            cntl = cntl + 1;
            sigsb5l(cntl, :, :) = sig(NgoodUnitsID, stimeind(i-1) - 1000 + 1:stimeind(i-1) + 5000);
    
        elseif b5o2(i)
    
            cnto2 = cnto2 + 1;
            sigsb5o2(cnto2, :, :) = sig(NgoodUnitsID, stimeind(i-1) - 1000 + 1:stimeind(i-1) + 5000);
    
        elseif b5o3(i)
    
            cnto3 = cnto3 + 1;
            sigsb5o3(cnto3, :, :) = sig(NgoodUnitsID, stimeind(i-2) - 1000 + 1:stimeind(i-2) + 5000);
    
        elseif b5o4(i)
    
            cnto4 = cnto4 + 1;
            sigsb5o4(cnto4, :, :) = sig(NgoodUnitsID, stimeind(i-3) - 1000 + 1:stimeind(i-3) + 5000);
    
        end
    
        if mod(i, 2000) == 1
            fprintf(num2str(i) + "-");
        end
    
    end

    disp("_");

    signalList = struct();
    signalList.area = areaName;

    signalList.sigsb2l = sigsb2l;
    signalList.sigsb2o2 = sigsb2o2;
    signalList.sigsb2o3 = sigsb2o3;
    signalList.sigsb2o4 = sigsb2o4;
    
    signalList.sigsb4l = sigsb4l;
    signalList.sigsb4o2 = sigsb4o2;
    signalList.sigsb4o3 = sigsb4o3;
    signalList.sigsb4o4 = sigsb4o4;
    
    signalList.sigsb5l = sigsb5l;
    signalList.sigsb5o2 = sigsb5o2;
    signalList.sigsb5o3 = sigsb5o3;
    signalList.sigsb5o4 = sigsb5o4;

end
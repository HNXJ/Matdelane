function OutPut = jGLOSignals(nwbPath, fileID, signalMode, time_pre, time_post, probeId, medFlag, chMed, twMed)

    disp(">Starting");
    
    try

        if medFlag

            fprintf("->2D Median filter is on with channel window radius = %d, time window radius %d\n", chMed, twMed);
            jFilter = "Median_chw_" + num2str(chMed) + "_twMed+" + num2str(twMed);
        
        else

            fprintf("->2D Median filter is off\n");
            jFilter = "off";

        end

    catch

        disp("->No proper 2D median filter parameters detected. Median filter is off.");

        medFlag = 0;
        chMed = 0;
        twMed = 0;
        jFilter = "off";

    end

    nwbFiles = {dir(nwbPath).name};
    nwbFiles = nwbFiles(endsWith(nwbFiles, ".nwb"));
    file_id = fileID;
    
    nwbFile = nwbPath + nwbFiles{file_id};
    nwb = nwbRead(nwbFile);
    ProbeNames = unique(nwb.general_extracellular_ephys_electrodes.vectordata.get('probe').data(:));
    ProbeArea = nwb.general_extracellular_ephys.get(ProbeNames{probeId+1}).location;

    disp("->Probe area is : " + string(ProbeArea));
    jGLOPassiveTrialCounts(nwb);
    a = nwb.acquisition.get("probe_" + num2str(probeId) + "_" + signalMode).electricalseries.get("probe_" + num2str(probeId) + "_" + signalMode + "_data");

    try
    
        e = nwb.intervals.get('passive_glo');
    
    catch

        disp("-> Warning! this session has no passive GLO mode; try using active function or check your file.");
        disp("> Exiting with NULL output");
        OutPut = [];
        return
    
    end
    
    eye = nwb.acquisition.get("eye_1_tracking").spatialseries.get("eye_1_tracking_data");
    pup = nwb.acquisition.get("pupil_1_tracking").timeseries.get("pupil_1_diameter_data").data(:);

    try
        rew = nwb.acquisition.get("reward_1_tracking").timeseries.get("reward_1_tracking_data").data(:);
    catch
        rew = pup * 0;
    end

    post_length = time_post;
    pre_length = time_pre;
    smoothness = 21;

    fprintf("\n->Processing trials - This function is ONLY for passive mode. For active, use <ActiveGLOSignals.m>\n");
    
    fs = a.starting_time_rate; % Sampling rate
    mua_signal = a.data(:, :); % Signals
    eye_signal = eye.data(:, :); % eye
    ch_count = size(mua_signal, 1); % Channel count
    
    try

        photodiode_signal = nwb.acquisition.get("photodiode_1_tracking").timeseries.get("photodiode_1_tracking_data").data(:);
   
    catch

        photodiode_signal = eye_signal(1, :) * 0;
        disp("This file does not have photodiode signal");

    end
    
    % Block 1 (HABEXP)
    lo_habexp = (e.vectordata.get('task_block_number').data(:) == 1) & (e.vectordata.get('stimulus_number').data(:) == 5) & e.vectordata.get('correct').data(:);
    
    % Block 2 (GLOEXP)
    go_gloexp = e.vectordata.get('go_gloexp').data(:) & e.vectordata.get('correct').data(:) & (e.vectordata.get('task_block_number').data(:) == 2);
    lo_gloexp = e.vectordata.get('lo_gloexp').data(:) & e.vectordata.get('correct').data(:) & (e.vectordata.get('task_block_number').data(:) == 2);
    
    % Block 3 (RNDCTL)
    go_rndctl = e.vectordata.get('go_rndctl').data(:) & e.vectordata.get('correct').data(:);
    lo_rndctl = e.vectordata.get('lo_rndctl').data(:) & e.vectordata.get('correct').data(:);
    igo_rndctl = e.vectordata.get('igo_rndctl').data(:) & e.vectordata.get('correct').data(:);
    ilo_rndctl = e.vectordata.get('ilo_rndctl').data(:) & e.vectordata.get('correct').data(:);
    
    % Block 4 (SEQCTL)
    go_seqctl = e.vectordata.get('go_seqctl').data(:) & e.vectordata.get('correct').data(:);
    igo_seqctl = e.vectordata.get('igo_seqctl').data(:) & e.vectordata.get('correct').data(:);

    correct_glo = e.vectordata.get('correct').data(:);
    % orientation_glo = e.vectordata.get('orientation').data(:);
    % presentation_glo = e.vectordata.get('presentation').data(:);

    t_start = e.start_time.data(:);
    % t_stop = e.stop_time.data(:);
    % t_stamp = (a.timestamps(:)-pre_length)*fs;
    
    r = ceil(post_length*fs);
    k = ceil(pre_length*fs);
    
    local_mua_hab = zeros(sum(lo_habexp), ch_count, r + k + 1); % 1 

    global_mua_ac = zeros(sum(go_gloexp), ch_count, r + k + 1); % 2
    local_mua_ac = zeros(sum(lo_gloexp), ch_count, r + k + 1); % 2

    g_ctl_mua_ac = zeros(sum(go_rndctl), ch_count, r + k + 1); % 3
    l_ctl_mua_ac = zeros(sum(lo_rndctl), ch_count, r + k + 1); % 3
    ig_ctl_mua_ac = zeros(sum(igo_rndctl), ch_count, r + k + 1); % 3
    il_ctl_mua_ac = zeros(sum(ilo_rndctl), ch_count, r + k + 1); % 3

    g_seq_mua_ac = zeros(sum(go_seqctl), ch_count, r + k + 1); % 4
    ig_seq_mua_ac = zeros(sum(igo_seqctl), ch_count, r + k + 1); % 4

    % global_mua_choice = zeros(1, ch_count, r + k + 1);
    % local_mua_choice = zeros(1, ch_count, r + k + 1);
    % global_mua_select = zeros(1, ch_count, r + k + 1);
    % local_mua_select = zeros(1, ch_count, r + k + 1);

    % global_c_mua_choice = zeros(1, ch_count, r + k + 1);
    % local_c_mua_choice = zeros(1, ch_count, r + k + 1);

    h_cnt = 0;
    g_cnt = 0;
    l_cnt = 0;

    gc_cnt = 0;
    lc_cnt = 0;
    igr_cnt = 0;
    ilr_cnt = 0;

    gs_cnt = 0;
    igs_cnt = 0;
    
    % gtemp_eye = zeros(1, k+r);
    % gtemp_eye2 = zeros(1, k+r);
    % gtemp_pup = zeros(1, k+r);
    gtemp_photodiode = zeros(1, k+r);
    
    g_reward = zeros(1, k+r);
    l_reward = zeros(1, k+r);
    
    % ltemp_eye = zeros(1, k+r);
    % ltemp_eye2 = zeros(1, k+r);
    % ltemp_pup = zeros(1, k+r);
    ltemp_photodiode = zeros(1, k+r);

    baseline_max = 1000;
    baseline_min = 500;

    % g_onset_choice = zeros(1, 1);
    % l_onset_choice = zeros(1, 1);
    % g_onset_select = zeros(1, 1);
    % l_onset_select = zeros(1, 1);
    % 
%         g_onset_reward = zeros(1, 1);
%         l_onset_reward = zeros(1, 1);
    % g_reaction_time = zeros(1, 1);
    % l_reaction_time = zeros(1, 1);

    interval_times = 1:size(t_start, 1);

    for i = interval_times
    
        start_index = floor(t_start(i) * fs);
        % stop_index = floor(t_stop(i) * fs);
        l = start_index;
    
        if correct_glo(i) == 0
    
            continue; % Skip incorrect trials
    
        end
    
        if lo_habexp(i)
    
            fprintf("-->B1:Local(C-habexp) %d\n", i);
            h_cnt = h_cnt + 1;
            
%                 gtemp_eye(h_cnt, :) = smooth(eye_signal(2, l-k+1:l+r), 10);
%                 gtemp_eye(h_cnt, :) = gtemp_eye(h_cnt, :) - min(gtemp_eye(h_cnt, 1:baseline_max));
%                 gtemp_eye(h_cnt, :) = gtemp_eye(h_cnt, :) / max(gtemp_eye(h_cnt, 1:baseline_max));
%     
            gtemp_photodiode(h_cnt, :) = smooth(photodiode_signal(l-k+1:l+r), 20);
%                 gtemp_photodiode(h_cnt, :) = gtemp_photodiode(h_cnt, :) - min(gtemp_photodiode(h_cnt, 1:baseline_max));
%                 gtemp_photodiode(h_cnt, :) = gtemp_photodiode(h_cnt, :) / max(gtemp_photodiode(h_cnt, 1:baseline_max));

            g_reward(h_cnt, :) = rew(l-k+1:l+r);
%                 g_onset_reward(h_cnt) = find(g_reward(h_cnt, :) > 2., 1);

%                 g_onset_choice(h_cnt) = find(gtemp_photodiode(h_cnt, 1000:end) < .5, 1) + 1000;
%                 g_onset_select(h_cnt) = find(gtemp_eye(h_cnt, 1200:end) > .5, 1) + 1200;
%                 g_reaction_time(h_cnt) = (g_onset_select(h_cnt) - g_onset_choice(h_cnt))/fs;

            if medFlag

                temp_sig = jMedianFilt2(mua_signal(:, l-k+1:l+r), chMed, twMed); 

            else

                temp_sig = mua_signal(:, l-k+1:l+r);

            end

%                 temp_sig_choice_aligned = mua_signal(:, 1+l+g_onset_choice(h_cnt)-1000:l+g_onset_choice(h_cnt)+r-1000);
%                 temp_sig_select_aligned = mua_signal(:, 1+l+g_onset_select(h_cnt)-1000:l+g_onset_select(h_cnt)+r-1000);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
%                     temp_sig_choice_aligned(ch, :) = smooth(temp_sig_choice_aligned(ch, :) - mean(temp_sig_choice_aligned(ch, 1:baseline_len)), smoothness);
%                     temp_sig_select_aligned(ch, :) = smooth(temp_sig_select_aligned(ch, :) - mean(temp_sig_select_aligned(ch, 1:baseline_len)), smoothness);
            end
    
            local_mua_hab(h_cnt, :, 1:r+k) = temp_sig;
%                 global_mua_choice(h_cnt, :, 1:r+k) = temp_sig_choice_aligned;
%                 global_mua_select(h_cnt, :, 1:r+k) = temp_sig_select_aligned;
    
        elseif go_gloexp(i)
    
            fprintf("-->B2:Global(C-gloexp) %d\n", i);
            g_cnt = g_cnt + 1;
            
%                 gtemp_eye(g_cnt, :) = smooth(eye_signal(2, l-k+1:l+r), 10);
%                 gtemp_eye(g_cnt, :) = gtemp_eye(g_cnt, :) - min(gtemp_eye(g_cnt, 1:baseline_max));
%                 gtemp_eye(g_cnt, :) = gtemp_eye(g_cnt, :) / max(gtemp_eye(g_cnt, 1:baseline_max));
%     
            gtemp_photodiode(g_cnt, :) = smooth(photodiode_signal(l-k+1:l+r), 20);
%                 gtemp_photodiode(g_cnt, :) = gtemp_photodiode(g_cnt, :) - min(gtemp_photodiode(g_cnt, 1:baseline_max));
%                 gtemp_photodiode(g_cnt, :) = gtemp_photodiode(g_cnt, :) / max(gtemp_photodiode(g_cnt, 1:baseline_max));

            g_reward(g_cnt, :) = rew(l-k+1:l+r);
%                 g_onset_reward(g_cnt) = find(g_reward(g_cnt, :) > 2., 1);

%                 g_onset_choice(g_cnt) = find(gtemp_photodiode(g_cnt, 1000:end) < .5, 1) + 1000;
%                 g_onset_select(g_cnt) = find(gtemp_eye(g_cnt, 1200:end) > .5, 1) + 1200;
%                 g_reaction_time(g_cnt) = (g_onset_select(g_cnt) - g_onset_choice(g_cnt))/fs;

            if medFlag 
                
                temp_sig = jMedianFilt2(mua_signal(:, l-k+1:l+r), chMed, twMed); 

            else

                temp_sig = mua_signal(:, l-k+1:l+r);

            end

%                 temp_sig_choice_aligned = mua_signal(:, 1+l+g_onset_choice(g_cnt)-1000:l+g_onset_choice(g_cnt)+r-1000);
%                 temp_sig_select_aligned = mua_signal(:, 1+l+g_onset_select(g_cnt)-1000:l+g_onset_select(g_cnt)+r-1000);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
%                     temp_sig_choice_aligned(ch, :) = smooth(temp_sig_choice_aligned(ch, :) - mean(temp_sig_choice_aligned(ch, 1:baseline_len)), smoothness);
%                     temp_sig_select_aligned(ch, :) = smooth(temp_sig_select_aligned(ch, :) - mean(temp_sig_select_aligned(ch, 1:baseline_len)), smoothness);
            end
    
            global_mua_ac(g_cnt, :, 1:r+k) = temp_sig;
%                 global_mua_choice(g_cnt, :, 1:r+k) = temp_sig_choice_aligned;
%                 global_mua_select(g_cnt, :, 1:r+k) = temp_sig_select_aligned;
    
        elseif lo_gloexp(i)
    
            fprintf("-->B2:Local(C-gloexp) %d\n", i);
            l_cnt = l_cnt + 1;
            
%                 ltemp_eye(l_cnt, :) = smooth(eye_signal(2, l-k+1:l+r), 10);
%                 ltemp_eye(l_cnt, :) = ltemp_eye(l_cnt, :) - min(ltemp_eye(l_cnt, 1:baseline_max));
%                 ltemp_eye(l_cnt, :) = ltemp_eye(l_cnt, :) / max(ltemp_eye(l_cnt, 1:baseline_max));

            ltemp_photodiode(l_cnt, :) = smooth(photodiode_signal(l-k+1:l+r), 20);
%                 ltemp_photodiode(l_cnt, :) = ltemp_photodiode(l_cnt, :) - min(ltemp_photodiode(l_cnt, 1:baseline_max));
%                 ltemp_photodiode(l_cnt, :) = ltemp_photodiode(l_cnt, :) / max(ltemp_photodiode(l_cnt, 1:baseline_max));
            
            l_reward(l_cnt, :) = rew(l-k+1:l+r);
%     
%                 try
%                     l_onset_reward(l_cnt) = find(l_reward(l_cnt, :) > 2., 1);
%                 catch
%                     l_onset_reward(l_cnt) = l_onset_reward(l_cnt-1);
%                     disp(i+.5);
%                 end

%                 l_onset_choice(l_cnt) = find(ltemp_photodiode(l_cnt, 1000:end) < .5, 1) + 1000;
%                 l_onset_select(l_cnt) = find(ltemp_eye(l_cnt, 1200:end) < .5, 1) + 1200;
%                 l_reaction_time(l_cnt) = (l_onset_select(l_cnt) - l_onset_choice(l_cnt))/fs;

            if medFlag 
                
                temp_sig = jMedianFilt2(mua_signal(:, l-k+1:l+r), chMed, twMed); 
             
            else

                temp_sig = mua_signal(:, l-k+1:l+r);
            
            end

%                 temp_sig_choice_aligned = mua_signal(:, 1+l+l_onset_choice(l_cnt)-1000:l+l_onset_choice(l_cnt)+r-1000);
%                 temp_sig_select_aligned = mua_signal(:, 1+l+l_onset_select(l_cnt)-1000:l+l_onset_select(l_cnt)+r-1000);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
%                     temp_sig_choice_aligned(ch, :) = smooth(temp_sig_choice_aligned(ch, :) - mean(temp_sig_choice_aligned(ch, 1:baseline_len)), smoothness);
%                     temp_sig_select_aligned(ch, :) = smooth(temp_sig_select_aligned(ch, :) - mean(temp_sig_select_aligned(ch, 1:baseline_len)), smoothness);
            end
    
            local_mua_ac(l_cnt, :, 1:r+k) = temp_sig;
%                 local_mua_choice(l_cnt, :, 1:r+k) = temp_sig_choice_aligned;
%                 local_mua_select(l_cnt, :, 1:r+k) = temp_sig_select_aligned;
    
        elseif go_rndctl(i)
        
            fprintf("-->B3:Global(C-rndctl) %d\n", i);            
            gc_cnt = gc_cnt + 1;
            temp_sig = mua_signal(:, l-k+1:l+r);
          
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
            end
    
            g_ctl_mua_ac(gc_cnt, :, 1:r+k) = temp_sig;
        
        elseif lo_rndctl(i)
        
            fprintf("-->B3:Local(C-rndctl) %d\n", i);
            lc_cnt = lc_cnt + 1;
            temp_sig = mua_signal(:, l-k+1:l+r);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
            end
    
            l_ctl_mua_ac(lc_cnt, :, 1:r+k) = temp_sig;
     
        elseif igo_rndctl(i)
    
            fprintf("-->B3:iGlobal(C-rndctl) %d\n", i);
            igr_cnt = igr_cnt + 1;
            
%                 gtemp_eye(igr_cnt, :) = smooth(eye_signal(2, l-k+1:l+r), 10);
%                 gtemp_eye(igr_cnt, :) = gtemp_eye(igr_cnt, :) - min(gtemp_eye(igr_cnt, 1:baseline_max));
%                 gtemp_eye(igr_cnt, :) = gtemp_eye(igr_cnt, :) / max(gtemp_eye(igr_cnt, 1:baseline_max));
%     
            gtemp_photodiode(igr_cnt, :) = smooth(photodiode_signal(l-k+1:l+r), 20);
%                 gtemp_photodiode(igr_cnt, :) = gtemp_photodiode(igr_cnt, :) - min(gtemp_photodiode(igr_cnt, 1:baseline_max));
%                 gtemp_photodiode(igr_cnt, :) = gtemp_photodiode(igr_cnt, :) / max(gtemp_photodiode(igr_cnt, 1:baseline_max));

            g_reward(igr_cnt, :) = rew(l-k+1:l+r);
%                 g_onset_reward(igr_cnt) = find(g_reward(igr_cnt, :) > 2., 1);

%                 g_onset_choice(igr_cnt) = find(gtemp_photodiode(igr_cnt, 1000:end) < .5, 1) + 1000;
%                 g_onset_select(igr_cnt) = find(gtemp_eye(igr_cnt, 1200:end) > .5, 1) + 1200;
%                 g_reaction_time(igr_cnt) = (g_onset_select(igr_cnt) - g_onset_choice(igr_cnt))/fs;

            if medFlag 
                
                temp_sig = jMedianFilt2(mua_signal(:, l-k+1:l+r), chMed, twMed); 
          
            else

                temp_sig = mua_signal(:, l-k+1:l+r);

            end

%                 temp_sig_choice_aligned = mua_signal(:, 1+l+g_onset_choice(igr_cnt)-1000:l+g_onset_choice(igr_cnt)+r-1000);
%                 temp_sig_select_aligned = mua_signal(:, 1+l+g_onset_select(igr_cnt)-1000:l+g_onset_select(igr_cnt)+r-1000);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
%                     temp_sig_choice_aligned(ch, :) = smooth(temp_sig_choice_aligned(ch, :) - mean(temp_sig_choice_aligned(ch, 1:baseline_len)), smoothness);
%                     temp_sig_select_aligned(ch, :) = smooth(temp_sig_select_aligned(ch, :) - mean(temp_sig_select_aligned(ch, 1:baseline_len)), smoothness);
            end
    
            ig_ctl_mua_ac(igr_cnt, :, 1:r+k) = temp_sig;
%                 global_mua_choice(igr_cnt, :, 1:r+k) = temp_sig_choice_aligned;
%                 global_mua_select(igr_cnt, :, 1:r+k) = temp_sig_select_aligned;

        elseif ilo_rndctl(i)
    
            fprintf("-->B3:iLocal(C-rndctl) %d\n", i);
            ilr_cnt = ilr_cnt + 1;
            
%                 gtemp_eye(ilr_cnt, :) = smooth(eye_signal(2, l-k+1:l+r), 10);
%                 gtemp_eye(ilr_cnt, :) = gtemp_eye(ilr_cnt, :) - min(gtemp_eye(ilr_cnt, 1:baseline_max));
%                 gtemp_eye(ilr_cnt, :) = gtemp_eye(ilr_cnt, :) / max(gtemp_eye(ilr_cnt, 1:baseline_max));
%     
            gtemp_photodiode(ilr_cnt, :) = smooth(photodiode_signal(l-k+1:l+r), 20);
%                 gtemp_photodiode(ilr_cnt, :) = gtemp_photodiode(ilr_cnt, :) - min(gtemp_photodiode(ilr_cnt, 1:baseline_max));
%                 gtemp_photodiode(ilr_cnt, :) = gtemp_photodiode(ilr_cnt, :) / max(gtemp_photodiode(ilr_cnt, 1:baseline_max));

            g_reward(ilr_cnt, :) = rew(l-k+1:l+r);
%                 g_onset_reward(ilr_cnt) = find(g_reward(ilr_cnt, :) > 2., 1);

%                 g_onset_choice(ilr_cnt) = find(gtemp_photodiode(ilr_cnt, 1000:end) < .5, 1) + 1000;
%                 g_onset_select(ilr_cnt) = find(gtemp_eye(ilr_cnt, 1200:end) > .5, 1) + 1200;
%                 g_reaction_time(ilr_cnt) = (g_onset_select(ilr_cnt) - g_onset_choice(ilr_cnt))/fs;

            if medFlag

                temp_sig = jMedianFilt2(mua_signal(:, l-k+1:l+r), chMed, twMed); 

            else

                temp_sig = mua_signal(:, l-k+1:l+r);

            end
%                 temp_sig_choice_aligned = mua_signal(:, 1+l+g_onset_choice(ilr_cnt)-1000:l+g_onset_choice(ilr_cnt)+r-1000);
%                 temp_sig_select_aligned = mua_signal(:, 1+l+g_onset_select(ilr_cnt)-1000:l+g_onset_select(ilr_cnt)+r-1000);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
%                     temp_sig_choice_aligned(ch, :) = smooth(temp_sig_choice_aligned(ch, :) - mean(temp_sig_choice_aligned(ch, 1:baseline_len)), smoothness);
%                     temp_sig_select_aligned(ch, :) = smooth(temp_sig_select_aligned(ch, :) - mean(temp_sig_select_aligned(ch, 1:baseline_len)), smoothness);
            end
    
            il_ctl_mua_ac(ilr_cnt, :, 1:r+k) = temp_sig;
%                 global_mua_choice(ilr_cnt, :, 1:r+k) = temp_sig_choice_aligned;
%                 global_mua_select(ilr_cnt, :, 1:r+k) = temp_sig_select_aligned;

        elseif go_seqctl(i)
    
            fprintf("-->B4:Global(C-seqctl) %d\n", i);
            gs_cnt = gs_cnt + 1;
            
%                 gtemp_eye(gs_cnt, :) = smooth(eye_signal(2, l-k+1:l+r), 10);
%                 gtemp_eye(gs_cnt, :) = gtemp_eye(gs_cnt, :) - min(gtemp_eye(gs_cnt, 1:baseline_max));
%                 gtemp_eye(gs_cnt, :) = gtemp_eye(gs_cnt, :) / max(gtemp_eye(gs_cnt, 1:baseline_max));
%     
            gtemp_photodiode(gs_cnt, :) = smooth(photodiode_signal(l-k+1:l+r), 20);
%                 gtemp_photodiode(gs_cnt, :) = gtemp_photodiode(gs_cnt, :) - min(gtemp_photodiode(gs_cnt, 1:baseline_max));
%                 gtemp_photodiode(gs_cnt, :) = gtemp_photodiode(gs_cnt, :) / max(gtemp_photodiode(gs_cnt, 1:baseline_max));

            g_reward(gs_cnt, :) = rew(l-k+1:l+r);
%                 g_onset_reward(gs_cnt) = find(g_reward(gs_cnt, :) > 2., 1);

%                 g_onset_choice(gs_cnt) = find(gtemp_photodiode(gs_cnt, 1000:end) < .5, 1) + 1000;
%                 g_onset_select(gs_cnt) = find(gtemp_eye(gs_cnt, 1200:end) > .5, 1) + 1200;
%                 g_reaction_time(gs_cnt) = (g_onset_select(gs_cnt) - g_onset_choice(gs_cnt))/fs;

            if medFlag 
                
                temp_sig = jMedianFilt2(mua_signal(:, l-k+1:l+r), chMed, twMed); 
            
            else

                temp_sig = mua_signal(:, l-k+1:l+r);

            end

%                 temp_sig_choice_aligned = mua_signal(:, 1+l+g_onset_choice(gs_cnt)-1000:l+g_onset_choice(gs_cnt)+r-1000);
%                 temp_sig_select_aligned = mua_signal(:, 1+l+g_onset_select(gs_cnt)-1000:l+g_onset_select(gs_cnt)+r-1000);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
%                     temp_sig_choice_aligned(ch, :) = smooth(temp_sig_choice_aligned(ch, :) - mean(temp_sig_choice_aligned(ch, 1:baseline_len)), smoothness);
%                     temp_sig_select_aligned(ch, :) = smooth(temp_sig_select_aligned(ch, :) - mean(temp_sig_select_aligned(ch, 1:baseline_len)), smoothness);
            end
    
            g_seq_mua_ac(gs_cnt, :, 1:r+k) = temp_sig;
%                 global_mua_choice(gs_cnt, :, 1:r+k) = temp_sig_choice_aligned;
%                 global_mua_select(gs_cnt, :, 1:r+k) = temp_sig_select_aligned;

        elseif igo_seqctl(i)
    
            fprintf("-->B4:iGlobal(C-seqctl) %d\n", i);
            igs_cnt = igs_cnt + 1;
            
%                 gtemp_eye(igs_cnt, :) = smooth(eye_signal(2, l-k+1:l+r), 10);
%                 gtemp_eye(igs_cnt, :) = gtemp_eye(igs_cnt, :) - min(gtemp_eye(igs_cnt, 1:baseline_max));
%                 gtemp_eye(igs_cnt, :) = gtemp_eye(igs_cnt, :) / max(gtemp_eye(igs_cnt, 1:baseline_max));
%     
            gtemp_photodiode(igs_cnt, :) = smooth(photodiode_signal(l-k+1:l+r), 20);
%                 gtemp_photodiode(igs_cnt, :) = gtemp_photodiode(igs_cnt, :) - min(gtemp_photodiode(igs_cnt, 1:baseline_max));
%                 gtemp_photodiode(igs_cnt, :) = gtemp_photodiode(igs_cnt, :) / max(gtemp_photodiode(igs_cnt, 1:baseline_max));

            g_reward(igs_cnt, :) = rew(l-k+1:l+r);
%                 g_onset_reward(igs_cnt) = find(g_reward(igs_cnt, :) > 2., 1);

%                 g_onset_choice(igs_cnt) = find(gtemp_photodiode(igs_cnt, 1000:end) < .5, 1) + 1000;
%                 g_onset_select(igs_cnt) = find(gtemp_eye(igs_cnt, 1200:end) > .5, 1) + 1200;
%                 g_reaction_time(igs_cnt) = (g_onset_select(igs_cnt) - g_onset_choice(igs_cnt))/fs;

            if medFlag 
                
                temp_sig = jMedianFilt2(mua_signal(:, l-k+1:l+r), chMed, twMed); 
            
            else

                temp_sig = mua_signal(:, l-k+1:l+r);

            end

%                 temp_sig_choice_aligned = mua_signal(:, 1+l+g_onset_choice(igs_cnt)-1000:l+g_onset_choice(igs_cnt)+r-1000);
%                 temp_sig_select_aligned = mua_signal(:, 1+l+g_onset_select(igs_cnt)-1000:l+g_onset_select(igs_cnt)+r-1000);
            
            for ch = 1:ch_count
                temp_sig(ch, :) = smooth(temp_sig(ch, :) - mean(temp_sig(ch, baseline_min:baseline_max)), smoothness);
%                     temp_sig_choice_aligned(ch, :) = smooth(temp_sig_choice_aligned(ch, :) - mean(temp_sig_choice_aligned(ch, 1:baseline_len)), smoothness);
%                     temp_sig_select_aligned(ch, :) = smooth(temp_sig_select_aligned(ch, :) - mean(temp_sig_select_aligned(ch, 1:baseline_len)), smoothness);
            end
    
            ig_seq_mua_ac(igs_cnt, :, 1:r+k) = temp_sig;
%                 global_mua_choice(igs_cnt, :, 1:r+k) = temp_sig_choice_aligned;
%                 global_mua_select(igs_cnt, :, 1:r+k) = temp_sig_select_aligned;

        end
    
    end
    
    disp("->Generating output struct ...");
    rk = r + k;

    OutPut = struct();

    OutPut.Fs = fs;
    OutPut.Session = string(nwb.identifier);
    OutPut.Area = string(ProbeArea);
    OutPut.Filter = jFilter;

    OutPut.B1L = local_mua_hab(:, :, 1:rk);   
    OutPut.B2L = local_mua_ac(:, :, 1:rk);
    OutPut.B2G = global_mua_ac(:, :, 1:rk);

    OutPut.B3L = l_ctl_mua_ac(:, :, 1:rk);
    OutPut.B3G = g_ctl_mua_ac(:, :, 1:rk);
    OutPut.B3iL = il_ctl_mua_ac(:, :, 1:rk);
    OutPut.B3iG = ig_ctl_mua_ac(:, :, 1:rk);
    
    OutPut.B4G = g_seq_mua_ac(:, :, 1:rk);
    OutPut.B4iG = ig_seq_mua_ac(:, :, 1:rk);

    OutPut.photodiode_g = gtemp_photodiode;
    OutPut.photodiode_l = ltemp_photodiode;
    OutPut.reward_g = g_reward;
    OutPut.reward_l = l_reward;

    disp(">Done." + string(nwb.identifier));
    
end

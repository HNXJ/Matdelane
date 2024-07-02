function chresp = jReceptiveFieldMap(nwb, p)

    ProbeNames = unique(nwb.general_extracellular_ephys_electrodes.vectordata.get('probe').data(:));
    numProbes = length(ProbeNames);
    identifier_ = nwb.identifier;
    ProbeAreas{p} = nwb.general_extracellular_ephys.get(ProbeNames{p}).location;
    
    fprintf("-> %d Probes detected in %s\n-->This function will process only probe no.%d (%s)", numProbes, identifier_, p, ProbeAreas{p}{1});
    
    start_times= nwb.intervals.get('rf_mapping_v2').start_time.data(:);
    rf_info=nwb.intervals.get('rf_mapping_v2');
    correct = rf_info.vectordata.get('correct').data(:);
    x_position = rf_info.vectordata.get('x_position').data(:);
    x_position_negative = rf_info.vectordata.get('x_position_negative').data(:);
    y_position = rf_info.vectordata.get('y_position').data(:);
    y_position_negative = rf_info.vectordata.get('y_position_negative').data(:); 
    sizes = rf_info.vectordata.get('size').data(:);
    % trials = rf_info.vectordata.get('trial_num').data(:);
    
    indx_correct_trls_from_all_events = find(~isnan(x_position) & correct);
    xpos = x_position;
    ypos = y_position;
    xnegs = find(x_position_negative);
    ynegs = find(y_position_negative);
    
    xpos(xnegs) = -xpos(xnegs);
    ypos(ynegs) = -ypos(ynegs);
    
    r = sizes;
    r = r(indx_correct_trls_from_all_events);
    xpos = xpos(indx_correct_trls_from_all_events);
    ypos = ypos(indx_correct_trls_from_all_events);
    
    conds = [xpos ypos r];
    uniconds = unique(conds,'rows');
    condindx = nan(1,size(conds,1));
    
    for cond = 1:length(uniconds)
        for c = 1:size(conds,1)
            if conds(c,:) == uniconds(cond,:)
                condindx(c)=cond;
            end
        end
    end
    
    mua_name = ['probe_',num2str(p-1),'_muae'];
    % PD = nwb.acquisition.get('photodiode_1_tracking').timeseries.get('photodiode_1_tracking_data');
    % PD_data = PD.data(:);
    mua = nwb.acquisition.get(mua_name).electricalseries.get([mua_name,'_data']);
    indx = jNearestIndex(mua.timestamps(:),start_times(indx_correct_trls_from_all_events));
    resp = jEpochData(mua.data(:,:),indx,[200,200]);
    resp = jBaselineCorrect(resp, 1:400);
    
    for chan = 1:size(resp, 1)
    
        for cond = 1:length(uniconds)
            chanresp(cond) = squeeze(mean(resp(chan,226+50:226+80,condindx==cond),[2,3]));
        end
    
        % chanresp(31) = 0;
        even_indx = 1:size(resp, 1);
        for cond = 1:length(uniconds)
            chanresp_grandavg(cond) = squeeze(mean(resp(even_indx,230:300,condindx==cond),[1,2,3]));
        end
        %chanresp_grandavg = chanresp_grandavg./max(chanresp_grandavg);
        chanresp_grandavg = zscore(chanresp_grandavg);
        % chanresp_grandavg(31) = 0;
        
        if mod(chan-1, 9) == 0
            figure("Position", [0 0 1800 1400]);
        end
    
        subplot(3, 3, mod(chan-1, 9)+1);
        bubblechart(uniconds(:,1),uniconds(:,2),uniconds(:,3),chanresp);
        bubblesize([4 27]);
        colormap(gca,"jet");
        % caxis([-3 3]);
        colorbar;
        title(num2str(chan));
        set(gcf,'position',[50 50 1750 1350]);
    
    end
    
    unix = uniconds(:,1);
    uniy = uniconds(:,2);
    
    xu = unix;
    yu = uniy;
    
    figure("Position", [0 0 1800 1200]);
    bubblechart(xu,yu,uniconds(:,3),chanresp_grandavg);
    colormap(gca,"jet");
    %caxis([0.7 1])
    colorbar;
    title("RF map for " + string(ProbeAreas{p}{1}));
    set(gcf,'position',[100 100 1300 1300]);

    chresp = chanresp;

end


function out_mat = jEpochData(in_mat, idxs, vrange)

    nd = ndims(in_mat);
    out_mat = [];
    for ii = 1 : numel(idxs)
        if size(in_mat,2)==1
            out_mat = cat(nd+1, out_mat, in_mat(idxs(ii)-vrange(1):idxs(ii)+vrange(2)));
        else
            out_mat = cat(nd+1, out_mat, in_mat(:,idxs(ii)-vrange(1):idxs(ii)+vrange(2)));
        end
    end

end


function dat = jBaselineCorrect(dat, bl_epoch)

    if ismatrix(dat)
        dat = dat - repmat(nanmean(dat(:,bl_epoch),2), 1, size(dat,2));
    elseif ndims(dat) == 3
        dat = dat - repmat(nanmean(dat(:,bl_epoch,:),2), 1, size(dat,2), 1);
    end

end


function I = jNearestIndex(aRef, aTest)
    
    if size(aTest,2) > 1
        aTest = aTest.';
    end
    
    if size(aRef,2) == 1
        aRef = aRef.';
    end
    
    % d = nan(numel(aTest), 1);
    % idx = nan(numel(aTest), 1);
    
    edges = [-Inf, mean([aRef(2:end); aRef(1:end-1)]), +Inf];
    I = discretize(aTest, edges);

end

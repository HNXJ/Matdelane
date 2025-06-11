%% SET PATH
    addpath('fieldtrip-20240119/'), ft_defaults
    % clear - recommended

%% PART 1: LOAD IN THE DATA
    %%
    [lfp, dataspike] = create_data(); % to change parameters adjust the function
    %% create a spike structure from the continuous structure
    spike = ft_checkdata(dataspike,'datatype', 'spike', 'feedback', 'yes');

%% PART 2: COMPUTE THE SPIKE TRIGGERED SPECTRA

    % 
    cfg = []; 
    cfg.foi = 10:100; 
    cfg.taper = 'hann'; 
    cfg.t_ftimwin = 5./cfg.foi;     
    sts = ft_spiketriggeredspectrum_convol(cfg,lfp,spike);

%% PART 3: COMPUTE MEASURE OF SPIKE-FIELD ASSOCIATION

    %% compute the measure of phase locking strength
    nNeurons = length(spike.label); 
    ppc_all = NaN(nNeurons, length(sts.freq)); 
    n_spikes = NaN(nNeurons, length(sts.freq)); 
    for iNeuron = 1:nNeurons

        if length(sts.time{iNeuron})<3, continue,end
        
        % the code to compute for a single neuron
        cfg = []; 
        cfg.spikechannel = iNeuron; 
        cfg.method = 'ppc1'; 
        stat_ppc = ft_spiketriggeredspectrum_stat(cfg, sts); % this function takes only one input at a time
        
        % collect the data in some array for all neurons
        ppc_all(iNeuron,:) = stat_ppc.ppc1; % collect the ppc values across units
        n_spikes(iNeuron,:) = stat_ppc.nspikes; % collect the number of spikes for each neuron
    end

%% PART 4: VISUALIZE INDIVIDUAL NEURONS AND COMPUTE THE PHASE HISTOGRAM

    % Note: this function still expects spike data to come in continuous
    % format, as one structure. 
    data_all = ft_appenddata([],lfp, dataspike); 
    sta_all = []; 
    for iNeuron = 1:nNeurons
        cfg = []; 
        cfg.spikechannel = dataspike.label{iNeuron};   
        cfg.channel = lfp.label{1};
        sta = ft_spiketriggeredaverage(cfg,data_all);
        sta_all(iNeuron,:) = sta.avg; 
    end   
    figure, imagesc(sta.time, 1:nNeurons, sta_all), colorbar
    xlabel('time')
    ylabel('neurons')

%% PART 5: VISUALIZE INDIVIDUAL NEURONS AND COMPUTE THE PHASE HISTOGRAM    
    % plot the angles, rayleigh test, ppc for one neuron
    for iNeuron = 1:10
        if length(sts.time{iNeuron})<3, continue,end
    
        f_ind = nearest(sts.freq, 40); 
        neuron_ind = iNeuron; 
        
        % compute the rayleigh test
        cfg = []; 
        cfg.spikechannel = neuron_ind; 
        cfg.method = 'ral';      
        stat_ral = ft_spiketriggeredspectrum_stat(cfg, sts); % this function takes only one input at a time
        
        % for illustration also show the plv
        cfg = []; 
        cfg.spikechannel = neuron_ind; 
        cfg.method = 'plv';      
        stat_plv = ft_spiketriggeredspectrum_stat(cfg, sts); % this function takes only one input at a time
               
        %
        figure, 
        subplot(2,2,1)
        % visualize the angles for this neuron
        angles = angle(sts.fourierspctrm{iNeuron}(:,:,f_ind));        
        rose(angles,20)        
        
        % visualize the sta
        subplot(2,2,2)
        plot(sta.time, sta_all(iNeuron,:),'k')
        xlabel('time')
        ylabel('sta')

        % visualize the ppc for this neuron
        subplot(2,2,3)
        plot(stat_ppc.freq, ppc_all(iNeuron,:),'k')
        xlabel('frequency')
        ylabel('ppc')
        
        subplot(2,2,4)
        plot(stat_plv.freq, stat_plv.plv,'k')
        xlabel('frequency')
        ylabel('plv')

        title(sprintf("Rayleigh test: %2.2f\n, num. spikes: %d", stat_ral.ral(f_ind),n_spikes(iNeuron,f_ind)))
    end

    
%% PART 6: COMPUTE GROUP STATISTICS FOR ALL THE NEURONS 

    % weighted average 
    [w_ppc,uw_ppc,thresh_ppc] = compute_avg_ppc(n_spikes,ppc_all);
    
    % plot the results 
    figure, 
    h(1)=plot(stat_ppc.freq, w_ppc,'g'); 
    hold on, 
    h(2)=plot(stat_ppc.freq, uw_ppc,'r'); 
    hold on, 
    h(3)=plot(stat_ppc.freq, thresh_ppc,'b');
    legend(h,'weighted', 'unweighted', 'threshold50')


%% PART 7: COMPUTE THE SPECTRUM IN A DIFFERENT WAY, TAKING THE FOURIER AROUND EACH SPIKE

    cfg = []; 
    cfg.timwin= [-0.1 0.1]; 
    cfg.foilim = [10 100]; 
    cfg.taper = 'hann'; 
    cfg.spikechannel = 1:length(spike.trial); 
    sts_fft = ft_spiketriggeredspectrum_fft(cfg, lfp, spike); 
    
    % 
    nNeurons = length(spike.label); 
    ppc_all = NaN(nNeurons, length(sts_fft.freq)); 
    n_spikes = NaN(nNeurons, length(sts_fft.freq)); 
    for iNeuron = 1:nNeurons

        if length(sts.time{iNeuron})<3, continue,end 
        
        % the code to compute for a single neuron
        cfg = []; 
        cfg.spikechannel = iNeuron; 
        cfg.method = 'ppc1'; 
        stat_ppc = ft_spiketriggeredspectrum_stat(cfg, sts_fft); % this function takes only one input at a time
        
        % collect the data in some array for all neurons
        ppc_all(iNeuron,:) = stat_ppc.ppc1; % collect the ppc values across units
        n_spikes(iNeuron,:) = stat_ppc.nspikes; % collect the number of spikes for each neuron
    end

    %
    [w_ppc,uw_ppc,thresh_ppc] = compute_avg_ppc(n_spikes,ppc_all);
    
    % plot the results 
    figure, 
    h(1)=plot(stat_ppc.freq, w_ppc,'g'); 
    hold on, 
    h(2)=plot(stat_ppc.freq, uw_ppc,'r'); 
    hold on, 
    h(3)=plot(stat_ppc.freq, thresh_ppc,'b');
    legend(h,'weighted', 'unweighted', 'threshold50')

%% PART 8: COMPUTE THE SPIKE FIELD COHERENCE
    cfg = []; 
    cfg.method = 'mtmfft'; 
    cfg.taper = 'dpss';
    cfg.foilim = [1 100];
    cfg.tapsmofrq = 8;
    cfg.output = 'fourier';
    ft = ft_freqanalysis(cfg, data_all);
    
    cfg = []; 
    cfg.method = 'coh';
    con = ft_connectivityanalysis(cfg, ft); 

    coh_spike_lfp = squeeze(con.cohspctrm(1,:,:));
    figure, plot(con.freq, nanmean(coh_spike_lfp,1)); 






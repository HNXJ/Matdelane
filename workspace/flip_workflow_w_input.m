clear all; close all

load('C:\Users\Sophyy\Box\Projects\passive_GLO\scripts\2023_9_MUA\info_with_manual_areas.mat');

areas = unique(info.area);
chan_sessions = cellfun(@(x) x(1:end-4), info.session, 'UniformOutput', false);
sessions = unique(chan_sessions);
probes = unique(strcat(chan_sessions,'_',info.probe,'_',info.putativearea));
probes = probes(1:45) %take out Le probes;
numProbes = length(probes);

save_fig_flag = true;
pre_epoch_flag = true;
pre_epoch_dir = 'E:\passiveGLO_epoched_by_probe';


% additional user vars if needed

data_dir = 'E:\passiveGLO\';
obj_dir = 'C:\Users\Sophyy\OneDrive - Vanderbilt\passive_GLO\graphs\250227_L4';

% initiate user input var
user_input.subsel_input = cell(1,numProbes);
user_input.ca_MUA = cell(1,numProbes);
user_input.ca_CSD = cell(1,numProbes);
user_input.manual_bad_chans = cell(1,numProbes);
user_input.split_probe = cell(1,numProbes);
user_input.throw_out_first = cell(1,numProbes);
user_input.throw_out_last = cell(1,numProbes);
user_input.L4 = cell(1,numProbes);
user_input.orientation = cell(1,numProbes);



% recording info from google sheet
ID = '1FSPJK3k2w1mdn-Efde5DYrnIINQTsq84deYzr3WZR20';
url_name = sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s', ID);
recording_info = webread(url_name);


for p = [11,13,15,19]%[38:45]%1:numProbes
    selected_probe = probes{p};

    clearvars nwb
    close all
    cd(data_dir)
    sess = selected_probe(1:strfind(selected_probe,'_probe')-1);
    probe = selected_probe(strfind(selected_probe,'probe'):strfind(selected_probe,'probe')+5);
    probe_num = double(selected_probe(strfind(selected_probe,'probe')+5)-64);
    probe_area = selected_probe(strfind(selected_probe,'probe')+7:end);

    curr_probe_indx = contains(info.session,sess) & contains(info.probe,probe);

    s_rec = find(contains(recording_info.Identifier,sess));
    spacings = str2num(recording_info.electrodeSpacing_um_{s_rec});
    tasks = recording_info.Experiment_Description{s_rec};

    if pre_epoch_flag
        cd(pre_epoch_dir)
        load([sess,'_',probe,'_',probe_area,'_epoched'])
        [sess,'_',probe,'_',probe_area]
    else
        [sess,'_',probe,'_',probe_area]
        nwb = nwbRead([data_dir,sess,'.nwb']);


        if contains(sess,'C31')
            PO_offset = 25;
        else
            PO_offset = 18;
        end

        %check if flash exists
        try
            start_times= nwb.intervals.get('flash').start_time.data(:);
            stims=nwb.intervals.get('flash').vectordata.get('stimulus_number').data(:);
            starts = stims == 2;
        catch
            start_times= nwb.intervals.get('passive_glo').start_time.data(:);
            presentation=nwb.intervals.get('passive_glo').vectordata.get('presentation').data(:);
            correct = nwb.intervals.get('passive_glo').vectordata.get('correct').data(:);
            starts = presentation == 1 & correct;
        end

        probeArea_nwb = nwb.general_extracellular_ephys.get(probe).location{:};
        if ~strcmp(probe_area,probeArea_nwb)
            disp("potential probe area discrepancy, check probe info")
        end


        lfp.ElecSeries = nwb.acquisition.get(['probe_',num2str(probe_num-1),'_lfp']).electricalseries.get(['probe_',num2str(probe_num-1),'_lfp','_data']);
        mua.ElecSeries = nwb.acquisition.get(['probe_',num2str(probe_num-1),'_muae']).electricalseries.get(['probe_',num2str(probe_num-1),'_muae','_data']);

        lfp.times = nearest_index(lfp.ElecSeries.timestamps(:),start_times(starts))+PO_offset;
        mua.times = nearest_index(mua.ElecSeries.timestamps(:),start_times(starts))+PO_offset;

        if p == 32
            lfp.epoched_flip = epoch_data(lfp.ElecSeries.data(:,:),lfp.times,[500,1]);
            lfp.epoched = epoch_data(lfp.ElecSeries.data(:,:),lfp.times,[1,500]);
            mua.epoched = epoch_data(mua.ElecSeries.data(:,:),mua.times,[100,500]);
            lfp.epoched_medfilted = medfilt1(lfp.epoched_flip,3,[],1);

            mua.baselined = baseline_correct(mua.epoched,[1:99]);

        else
            lfp.epoched = epoch_data(lfp.ElecSeries.data(:,:),lfp.times,[1,500]);
            mua.epoched = epoch_data(mua.ElecSeries.data(:,:),mua.times,[100,500]);
            lfp.epoched_medfilted = medfilt1(lfp.epoched,3,[],1);
            mua.baselined = baseline_correct(mua.epoched,[1:99]);
        end

    end
    %select every other
    chansel.Nchan = size(lfp.epoched,1);
    chansel.all_idx = [1:chansel.Nchan];
    chansel.zig_idx= [1:2:chansel.Nchan];
    chansel.zag_idx= chansel.zig_idx+1;

    sig=0.4; %cortical conductivity
    dis = spacings(probe_num); %inter-site distance in micrometers
    sep = 2; %separation between channels
    chansel.label_sep = 2;
    if dis<50
        sep = 5;
        chansel.label_sep = 3;
    end

    close all

    %setting up yticks and yticklabels
    chansel.all_chan_num_str= num2str([1:chansel.Nchan]');
    chansel.all_labels = num2str([0:dis:dis*(chansel.Nchan-1)]');
    chansel.all_labels= [chansel.all_chan_num_str, repmat('|',chansel.Nchan,1) ,chansel.all_labels];
    chansel.zig_labels= chansel.all_labels(chansel.zig_idx,:);
    chansel.zag_labels= chansel.all_labels(chansel.zag_idx,:);

    chansel.ytix_main = chansel.all_idx(1:chansel.label_sep:end);
    chansel.ytixlb_main = chansel.all_labels(1:chansel.label_sep:end,:);
    chansel.ytix_zig = 1:chansel.label_sep:length(chansel.zig_idx);
    chansel.ytixlb_zig = chansel.zig_labels(1:chansel.label_sep:end,:);
    chansel.ytix_zag = 1:chansel.label_sep:length(chansel.zag_idx);
    chansel.ytixlb_zag = chansel.zag_labels(1:chansel.label_sep:end,:);



    %results


    flipObj.all = vFLIP2(lfp.epoched_medfilted(:,:,:), 'DataType', 'raw_cut','fsample', 1000, 'intdist', dis/1000);
    flipObj.set1 = vFLIP2(lfp.epoched_medfilted(chansel.zig_idx,:,:), 'DataType', 'raw_cut','fsample', 1000, 'intdist', dis/1000);
    flipObj.set2 = vFLIP2(lfp.epoched_medfilted(chansel.zag_idx,:,:), 'DataType', 'raw_cut','fsample', 1000, 'intdist', dis/1000);

    flipObj.all.Results
    flipObj.all.plot_result()
    yticks(chansel.ytix_main)
    yticklabels(chansel.ytixlb_main)
    ylabel('chanNum | depth (um)')
    sgtitle([sess,' ',probe,' ',probe_area,' whole probe flip'])
    set(gcf,'position',[0,750,950,565])

    flipObj.set1.Results
    flipObj.set1.plot_result()
    yticks(chansel.ytix_zig)
    yticklabels(chansel.ytixlb_zig)
    ylabel('chanNum | depth (um)')
    sgtitle([sess,' ',probe,' ',probe_area,' zig flip'])
    set(gcf,'position',[900,750,950,565])

    flipObj.set2.Results
    flipObj.set2.plot_result()
    yticks(chansel.ytix_zag)
    yticklabels(chansel.ytixlb_zag)
    ylabel('chanNum | depth (um)')
    sgtitle([sess,' ',probe,' ',probe_area,' zag flip'])
    set(gcf,'position',[900,100,950,565])


    user_input.subsel_input{p} = input('all (1), zig (2), or zag(3)? ');

    switch user_input.subsel_input{p}
        case 1
            chansel.zigzag_sel = chansel.all_idx;
            chansel.zsel_ytix = chansel.ytix_main;
            chansel.zsel_ytixlb = chansel.ytixlb_main;
        case 2
            chansel.zigzag_sel = chansel.zig_idx;
            chansel.zsel_ytix = chansel.ytix_zig;
            chansel.zsel_ytixlb = chansel.ytixlb_zig;
        case 3
            chansel.zigzag_sel = chansel.zag_idx;
            chansel.zsel_ytix = chansel.ytix_zag;
            chansel.zsel_ytixlb = chansel.ytixlb_zag;
    end

    user_input.manual_bad_chans{p} = input('any channels to interp? input (native) channels nums in array form, or nan: ');

    manual_good_chansel = ones(1,chansel.Nchan);
    if ~isnan(user_input.manual_bad_chans{p})
        manual_good_chansel(user_input.manual_bad_chans{p})=0;
    end
    chansel.zseled_good_chans = find(manual_good_chansel(chansel.zigzag_sel));


    lfp.selected = lfp.epoched_medfilted(chansel.zigzag_sel,:,:);
    lfp.selected_unfiltered = lfp.epoched(chansel.zigzag_sel,:,:);

    NTrl = size(lfp.selected,3);
    NSampl = size(lfp.selected,2);
    NCh_zsel = size(lfp.selected,1);

    sample_points = [1:NCh_zsel];
    sample_points = sample_points(chansel.zseled_good_chans);

    lfp.interped = lfp.selected;
    lfp.interped_unfiltered = lfp.selected;

    if ~isnan(user_input.manual_bad_chans{p})
        %spline interpolation
        for trl = 1:NTrl
            lfp.interped_unfiltered(:,:,trl) = interp1(sample_points,lfp.selected_unfiltered(chansel.zseled_good_chans,:,trl),[1:NCh_zsel],'spline');
            lfp.interped(:,:,trl) = interp1(sample_points,lfp.selected(chansel.zseled_good_chans,:,trl),[1:NCh_zsel],'spline');
        end
    end

    %CSD calc


    CSD1=nan(NCh_zsel,NSampl,NTrl);

    for ch = sep+1:NCh_zsel-sep
        CSD1(ch,:,:) =-sig*(lfp.interped_unfiltered(ch-sep,:,:)-(2.*lfp.interped_unfiltered(ch,:,:))+lfp.interped_unfiltered(ch+sep,:,:))/(((dis*0.001).*sep).^2);
    end


    L4_good = 0;
    ct = 1;
    while L4_good ==0

        % if ct>1
        close all

        flipObj.all.Results
        flipObj.all.plot_result()
        yticks(chansel.ytix_main)
        yticklabels(chansel.ytixlb_main)
        ylabel('chanNum | depth (um)')
        sgtitle([sess,' ',probe,' ',probe_area,' whole probe flip'])
        set(gcf,'position',[0,750,950,565])

        flipObj.zsel = vFLIP2(lfp.interped, 'DataType', 'raw_cut','fsample', 1000, 'intdist', dis/1000);
        flipObj.zsel.plot_result()
        yticks(chansel.zsel_ytix)
        yticklabels(chansel.zsel_ytixlb)
        ylabel('chanNum | depth (um)')
        sgtitle([sess,' ',probe,' ',probe_area,' selected+smoothed flip'])
        set(gcf,'position',[900,750,950,565])

        user_input.manual_bands_flag{p} = input('manual band selection? (1 for yes, 0 for no) ');

        if user_input.manual_bands_flag{p}
            user_input.manual_bands_good{p} = 0;
            while   user_input.manual_bands_good{p} == 0

                user_input.zsel_manual_hi_band{p} = input('high band: ');
                user_input.zsel_manual_lo_band{p} = input('low band: ');

                close all

                flipObj.all = vFLIP2(lfp.epoched_medfilted(:,:,:), 'DataType', 'raw_cut','fsample', 1000, ...
                    'intdist', dis/1000,'manual_hi_band',user_input.zsel_manual_hi_band{p},'manual_lo_band',user_input.zsel_manual_lo_band{p});
                flipObj.all.Results
                flipObj.all.plot_result()
                yticks(chansel.ytix_main)
                yticklabels(chansel.ytixlb_main)
                ylabel('chanNum | depth (um)')
                sgtitle([sess,' ',probe,' ',probe_area,' whole probe flip'])
                set(gcf,'position',[0,750,950,565])

                flipObj.zsel = vFLIP2(lfp.interped(:,:,:), 'DataType', 'raw_cut','fsample', 1000, ...
                    'intdist', dis/1000,'manual_hi_band',user_input.zsel_manual_hi_band{p},'manual_lo_band',user_input.zsel_manual_lo_band{p});
                flipObj.zsel.plot_result()
                yticks(chansel.zsel_ytix)
                yticklabels(chansel.zsel_ytixlb)
                ylabel('chanNum | depth (um)')
                sgtitle([sess,' ',probe,' ',probe_area,' selected+smoothed flip'])
                set(gcf,'position',[900,750,950,565])

                user_input.manual_bands_good{p} = input('manual band good? ');

            end

        end
        % end

        figure(5)
        set(gcf,'position',[1900,500,650,800])

        tiledlayout(2,1)

        nexttile(1)
        imagesc(squeeze(mean(mua.baselined,3)))
        xline(100,'r','LineWidth',2)
        yticks(chansel.ytix_main)
        yticklabels(chansel.ytixlb_main)
        colorbar
        if ct>1
            caxis(user_input.ca_MUA{p})
        end
        title('mua baselined')
        ylabel('chanNum | depth (um)')

        nexttile(2)
        imagesc(squeeze(mean(CSD1,3)))
        yticks(chansel.zsel_ytix)
        yticklabels(chansel.zsel_ytixlb)
        xlim([1 150])
        colorbar
        if ct>1
            caxis(user_input.ca_CSD{p})
        end
        title(['CSD | neighbor = ',num2str(sep)])
        ylabel('chanNum | depth (um)')

        ca_good = input('caxis good? ');

        while ca_good == 0

            user_input.ca_MUA{p} = input('caxis for MUA? ');
            nexttile(1)
            caxis(user_input.ca_MUA{p})
            user_input.ca_CSD{p} = input('caxis for CSD? ');
            nexttile(2)
            caxis(user_input.ca_CSD{p})
            ca_good = input('caxis good? ');
        end



        RF_start_chans = [];

        % get RF channel references
        for a = [2:5,7:12,15]
            curr_area_indx = contains(info.area(curr_probe_indx), areas{a});

            if sum(curr_area_indx) > 3
                curr_area_indx = find(curr_area_indx);
                startChan = curr_area_indx(1);
                RF_start_chans = [RF_start_chans,startChan];
            end
        end

        if ct == 1
            user_input.split_probe{p} = input(['split probe? 0 for no, (native) channel number for first chan of second area. RF suggests ',num2str(RF_start_chans),' : ']);
            user_input.throw_out_first{p} = input('is top split bad? ');
            user_input.throw_out_last{p} = input('is bottom split bad? ');
        end

        if user_input.split_probe{p} == 0
            Nsplits = 1;
        else
            Nsplits = length(user_input.split_probe{p})+1;
        end

        split_chans = cell(1,Nsplits);

        for s = 1:Nsplits-1
            [~, chansel.zsel_split{s}] = min(abs(chansel.zigzag_sel - user_input.split_probe{p}(s)));
        end

        switch Nsplits
            case 1
                split_chans{1} = [1:NCh_zsel];

            case 2
                split_chans{1} = [1:chansel.zsel_split{1}-1];
                split_chans{2} = [chansel.zsel_split{1}:NCh_zsel];
            case 3
                split_chans{1} = [1:chansel.zsel_split{1}-1];
                split_chans{2} = [chansel.zsel_split{1}:chansel.zsel_split{2}-1];
                split_chans{3} = [chansel.zsel_split{2}:NCh_zsel];
        end

        % close all
        chansel.zsel_all_ytixlb = chansel.all_labels(chansel.zigzag_sel,:);

        for s = 1:Nsplits
            if length(split_chans{s})<20
                chansel.split_label_sep{s} = 1;
            elseif length(split_chans{s})<50
                chansel.split_label_sep{s} = 2;
            else
                chansel.split_label_sep{s} = 3;
            end


            chansel.split_chans{s} = chansel.zigzag_sel(split_chans{s});

            chansel.split_chans_ytixlb{s} = chansel.zsel_all_ytixlb(split_chans{s},:);

            if user_input.throw_out_first{p} && s == 1
                continue
            end

            if user_input.throw_out_last{p} && s == Nsplits
                continue
            end

            figure(s+5)
            flipObj.final{s} = vFLIP2(lfp.interped(split_chans{s},:,:), 'DataType', 'raw_cut','fsample', 1000, 'intdist', dis/1000);
            flipObj.final{s}.Results
            flipObj.final{s}.plot_result()
            ylabel('chanNum | depth (um)')
            yticks([1:chansel.split_label_sep{s}:length(split_chans{s})])
            yticklabels(chansel.split_chans_ytixlb{s}(1:chansel.split_label_sep{s}:end,:))

            set(gcf,'position',[0,100+20*(s-1),950,565])

            user_input.manual_bands_flag_splits{p}{s} = input('manual band selection for this split? (1 for yes, 0 for no) ');

            if user_input.manual_bands_flag_splits{p}{s}
                user_input.manual_bands_flag_splits_good{p}{s} = 0;
                while   user_input.manual_bands_flag_splits_good{p}{s} == 0

                    user_input.zsel_manual_hi_band_splits{p}{s} = input('high band: ');
                    user_input.zsel_manual_lo_band_splits{p}{s} = input('low band: ');


                    flipObj.final{s} = vFLIP2(lfp.interped(split_chans{s},:,:), 'DataType', 'raw_cut','fsample', 1000, ...
                        'intdist', dis/1000,'manual_hi_band',user_input.zsel_manual_hi_band_splits{p}{s},'manual_lo_band',user_input.zsel_manual_lo_band_splits{p}{s});
                    flipObj.final{s}.Results
                    flipObj.final{s}.plot_result()
                    ylabel('chanNum | depth (um)')
                    yticks([1:chansel.split_label_sep{s}:length(split_chans{s})])
                    yticklabels(chansel.split_chans_ytixlb{s}(1:chansel.split_label_sep{s}:end,:))

                    user_input.manual_bands_flag_splits_good{p}{s} = input('manual band for this split good? ');

                end


            end

                            if ~isnan(flipObj.final{s}.Results.crossoverchannel)
                    L4_flip = chansel.split_chans_ytixlb{s}(flipObj.final{s}.Results.crossoverchannel);
                    display(['native flip shows L4 as ', num2str(L4_flip)])
                else
                    display(['native flip failed'])
                end

                if save_fig_flag
                    save_fig([sess,'_',probe,'_',probe_area,'_flip',num2str(s)],obj_dir)
                end
        end

        user_input.L4{p} = input('L4? ');
        figure(5)
        nexttile(1)
        yline(user_input.L4{p},'LineWidth',1,'Color','r')
        nexttile(2)
        clearvars zsel_L4
        for l = 1:length(user_input.L4{p})
            [~, zsel_L4(l)] = min(abs(chansel.zigzag_sel - user_input.L4{p}(l)))
        end
        yline(zsel_L4,'LineWidth',1,'Color','r')

        figure(1)
        yline(user_input.L4{p},'LineWidth',1,'Color','r')


        figure(2)
        yline(zsel_L4,'LineWidth',1,'Color','r')

        L4_good = input('L4 good? ');
        ct = ct+1;
    end

    user_input.orientation{s} = input('orientation for area(s)? 1 = sup-deep; 2 = deep-sup ');


    if save_fig_flag
        figure(5)
        save_fig([sess,'_',probe,'_',probe_area,'_MUA_CSD'],obj_dir)
        figure(1)
        save_fig([sess,'_',probe,'_',probe_area,'_raw_flip'],obj_dir)
        figure(2)
        save_fig([sess,'_',probe,'_',probe_area,'_zsel_flip'],obj_dir)
    end

    % pause


end




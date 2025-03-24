%%
%data dir
cd('D:\PropofolOddball')

state = 2; 
switch state
    case 1 
        %awake
        filestoread = dir('*-preOddball-*');
        load('052422_alphabeta_bins_mat_all_awake.mat')
    case 2
        %anes
        filestoread = dir('*-Anesthesia-*');
        load('052422_alphabeta_bins_mat_all_anesthesia.mat')
end

%load('D:\PropofolOddball\lfpmuaTrial\original\051922_chansel_for_sw.mat')
num_bins = 8
%housekeeping vars
trl_phase_mat{1,size(filestoread,1)+1} = nan;
tone_phase_mat{1,size(filestoread,1)} = nan;
tone_times_mat{1,size(filestoread,1)} = nan;
step_size = 2*pi/num_bins;
ta1=[20,28,30,4,28,22,15,15]; %defining time bounds of anes states (in mins)
ta2=[96,105,80,41,77,67,42,33];
%%
%pipeline

%get mua trials
clear lfp trialInfo
tic
cd('D:\PropofolOddball')
switch state
    case 1 
    %awake
        filestoread = dir('*-preOddball-*');
    case 2
    %anes
        filestoread = dir('*-Anesthesia-*');
end

num_bins=8
areas = {'Tpt','FEF'}%,'7','PFC'}
mua_by_phase_sess = nan(8,num_bins,2);%sess * phase * area
for sess = 1:size(filestoread,1)
    load(filestoread(sess).name,'amua','electrodeInfo');
    switch state
        case 1 
            %awake
            bins = bins_mat_all_ab_awake{sess};
        case 2
            %anes
            t1_ms = ta1(sess)*60*1000;
            t2_ms = ta2(sess)*60*1000;
            amua = amua(t1_ms:t2_ms,:);
            bins = bins_mat_all_ab_anesthesia{sess};
    end

    for a = 1:2
        if a==1
            area_indx=find(contains(electrodeInfo.area,'Tpt')|contains(electrodeInfo.area,'CPB'));
        else
            area_indx=find(contains(electrodeInfo.area,'FEF'));
        end
        area_bins = bins(:,area_indx);
        area_mua = amua(:,area_indx);
        
        for b = 1:num_bins
            mua_by_phase_sess(sess,b,a)=squeeze(mean(area_mua(area_bins==b),[1,2]));
        end
    end
end

mua_by_phase = squeeze(mean(mua_by_phase_sess,1,'omitnan'));

toc
%{
switch state
    case 1 
        %awake
        cd('D:\PropofolOddball\alphabeta_mua\awake')
    case 2
        %anes
        cd('D:\PropofolOddball\alphabeta_mua\anes')
        save(['090122_mua_avg_by_ab_phase.mat'],'mua_by_phase','-v7.3')
end
%%

%stats and graphing bins
a=1
c = jet(num_bins);
step_size = 2*pi/num_bins;

fit1 = fit((1:num_bins)',mua_by_phase(:,a),'fourier1','Upper',[Inf Inf Inf step_size],'Lower',[-Inf -Inf -Inf step_size]);
resids = nan(1,num_bins);
sst_mat = nan(1,num_bins);
for i=1:num_bins
    resids(i) = (mua_by_phase(i,a)-fit1(i))^2;
    sst_mat(i) = (mua_by_phase(i,a)-mean(mua_by_phase(:,a)))^2;
end
r2 = 1-(sum(resids)/sum(sst_mat));

bin_names{num_bins}=nan;
for i = 1:num_bins
    lower_bound = -pi+step_size*(i-1);
    higher_bound = -pi+step_size*i;
    bin_names{i}=[num2str(lower_bound),' to ',num2str(higher_bound)]
end
figure
for i = 1:num_bins
    bar(i,mua_by_phase(i,a),'FaceColor',c(i,:))
    hold on
end
plot(fit1)
%ylim([0.0045,0.0049]) %anes tpt
%ylim([0.0051,0.00516]) %anes fef
%ylim([0.0055,0.0058])%awake tpt
ylim([0.00553,0.00555])%awake tpt


hold off


title(['Average MUA by alpha/beta phase bin in ',areas{a},'| r2 = ',num2str(r2)])
fn = ['090122_mua_avg_by_ab_phase_',areas{a}]
switch state
    case 1 
        %awake
        cd('C:\Users\Sophyy\Box\Projects\Oddball_2021\graphs\090122_tones_by_ab_phase_awake')
    case 2
        %anes
        cd('C:\Users\Sophyy\Box\Projects\Oddball_2021\graphs\090122_tones_by_ab_phase_anes')
end

print([fn,'.eps'],'-depsc')
%save([fn,'.mat'])
saveas(gcf,[fn,'.jpg'])
saveas(gcf,[fn,'.fig'])
%}
%%
%stupid code to save fit stats

a=2
c = jet(num_bins);
step_size = 2*pi/num_bins;

fit1= fit((1:num_bins)',mua_by_phase(:,a),'fourier1','Upper',[Inf Inf Inf step_size],'Lower',[-Inf -Inf -Inf step_size]);
fit_mat{state,a}= fit1;
resids = nan(1,num_bins);
sst_mat = nan(1,num_bins);
for i=1:num_bins
    resids(i) = (mua_by_phase(i,a)-fit1(i))^2;
    sst_mat(i) = (mua_by_phase(i,a)-mean(mua_by_phase(:,a)))^2;
end
r2_mat{state,a} = 1-(sum(resids)/sum(sst_mat));

cd('C:\Users\Sophyy\Box\Projects\Oddball_2021\process_data')

save('mua_by_ab_phase_fitstats.mat','fit_mat','r2_mat')

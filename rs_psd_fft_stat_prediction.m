% Compute statistical contrast of survivors and non-survivors as well as evluation of prediction based on training and test sets
%  - requires input of subjects from rs_psd_subjects_for_study.m
%%
addpath 'D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\scripts\finalized_resuscitation\dependencies'
addpath 'D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\scripts\finalized_resuscitation\dependencies\altmany-export_fig-9ac0917'
addpath 'D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\scripts\finalized_resuscitation\dependencies\ft_boundedline\'
addpath 'D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\scripts\finalized_resuscitation\dependencies\subaxis'
addpath 'D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\scripts\finalized_resuscitation\dependencies\cbrewer'
addpath 'C:\Users\tkusterm\Documents\FieldTrip\fieldtrip-20180702'
ft_defaults
addpath 'C:\Users\tkusterm\Documents\FieldTrip\fieldtrip-20180702\external\brewermap'
%% backup of figures before running the script
keep_old('D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\Figures\Spectra','D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\Figures\Spectra\old\')
%% Select groups based on day of interest
% path that will lead to folder where patient data are stored
base_path = 'D:\Lausanne_analyses\RestingEEG\resting_data\';

% Day 1
training_pos = training_pos_d1;
training_neg = training_neg_d1;
test_pos = test_pos_d1;
test_neg = test_neg_d1;
% Day 2
% training_pos = training_pos_d2;
% training_neg = training_neg_d2;
% test_pos = test_pos_d2;
% test_neg = test_neg_d2;
keep training_pos training_neg test_pos test_neg base_path
%% load freq spectrum by outcome group
file_name = 'data_freq_mtm.mat';
% TRAINING
% pos outcome
for subj_iter = 1:numel(training_pos)
    file_to_load=[base_path,training_pos{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    % select frequencies and channels
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird');
    cfg.channel= {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_pos{subj_iter} = frq;
    % normalize power spectra for each subject
    freq_all_norm_pos{subj_iter} = freq_all_pos{subj_iter};
    freq_all_norm_pos{subj_iter}.powspctrm = bsxfun(@rdivide, freq_all_pos{subj_iter}.powspctrm, sum(freq_all_pos{subj_iter}.powspctrm,2));
end;
% neg outcome
for subj_iter = 1:numel(training_neg)
    file_to_load=[base_path,training_neg{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    % select frequencies and channels
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird'); 
    cfg.channel = {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_neg{subj_iter} = frq;
    % normalize power spectra for each subject
    freq_all_norm_neg{subj_iter} = freq_all_neg{subj_iter};
    freq_all_norm_neg{subj_iter}.powspctrm = bsxfun(@rdivide, freq_all_neg{subj_iter}.powspctrm, sum(freq_all_neg{subj_iter}.powspctrm,2));
end;

% TEST
% pos outcome
for subj_iter = 1:numel(test_pos)
    file_to_load=[base_path,test_pos{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    % select frequencies and channels
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird'); 
    cfg.channel= {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_pos_test{subj_iter} = frq;
    % normalize power spectra for each subject
    freq_all_norm_pos_test{subj_iter} = freq_all_pos_test{subj_iter};
    freq_all_norm_pos_test{subj_iter}.powspctrm = bsxfun(@rdivide, freq_all_pos_test{subj_iter}.powspctrm, sum(freq_all_pos_test{subj_iter}.powspctrm,2));
end;
% neg outcome
for subj_iter = 1:numel(test_neg)
    file_to_load=[base_path,test_neg{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    % select frequencies and channels
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird'); 
    cfg.channel = {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_neg_test{subj_iter} = frq;
    % normalize power spectra for each subject
    freq_all_norm_neg_test{subj_iter} = freq_all_neg_test{subj_iter};
    freq_all_norm_neg_test{subj_iter}.powspctrm = bsxfun(@rdivide, freq_all_neg_test{subj_iter}.powspctrm, sum(freq_all_neg_test{subj_iter}.powspctrm,2));
end;
%% grand average the data
cfg=[];
cfg.keepindividual = 'yes';

% training
freq_pos_ga=ft_freqgrandaverage(cfg,freq_all_pos{:});
freq_norm_pos_ga=ft_freqgrandaverage(cfg,freq_all_norm_pos{:});
freq_neg_ga=ft_freqgrandaverage(cfg,freq_all_neg{:});
freq_norm_neg_ga=ft_freqgrandaverage(cfg,freq_all_norm_neg{:});

% test
freq_pos_test_ga=ft_freqgrandaverage(cfg,freq_all_pos_test{:});
freq_norm_pos_test_ga=ft_freqgrandaverage(cfg,freq_all_norm_pos_test{:});
freq_neg_test_ga=ft_freqgrandaverage(cfg,freq_all_neg_test{:});
freq_norm_neg_test_ga=ft_freqgrandaverage(cfg,freq_all_norm_neg_test{:});
%% concatenate across groups to plot individuals
freq_all = {freq_all_pos{:},freq_all_pos_test{:},freq_all_neg{:},freq_all_neg_test{:}};
freq_all_norm = {freq_all_norm_pos{:},freq_all_norm_pos_test{:},freq_all_norm_neg{:},freq_all_norm_neg_test{:}};
subjects = {training_pos{:},test_pos{:},training_neg{:},test_neg{:}};
%% Plot all individiual subjects - spectrum only - normalized
fig1=figure;
plot_iter = 1;
ft_hastoolbox('brewermap', 1)
for i = 1:numel(freq_all_norm)
    subplot(5,ceil(numel(freq_all_norm)/5),plot_iter)
    plot_iter = plot_iter+1;
    cfg=[];
    cfg.xlim = [5 15];
    cfg.ylim = [0 0.05];
    cfg.layout='standard_1020.elc';
    elecs=p_layout('ladybird'); 
    cfg.channel = {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}};
    ft_singleplotER(cfg,freq_all_norm{i})
%     title([subjects{i}(1:4),', #: ',sprintf('%i',i)]);
    title(subjects{i}(1:4));
    if i ~= 1
        set(gca,'ytick',[])
    end
end;
mtit('Normalized spectra')
fig1.Color=[1 1 1];
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% export_fig('Z:\\USERS_DATA\\tkusterm\\resting_state\\Figures\\nothresh\\spectral\\prediction\\d1_individual_spectra_10Hz.png','-painters','-r300')
% export_fig('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_individual_spectra_10Hz.png','-painters','-r300')
%% plot group-wise normalized spectrum and topographies using ft_boundedline
foi=[7 14]; % frequencies to be plotted on topographies
ft_hastoolbox('brewermap',1) % color palette used
fig1=figure;
% plot the spectra in a lineplot ( plot 1/2)
subplot(1,4,[1 2])
cfg=[];
cfg.frequency = [1 40];
% cfg.ylim=[0 0.01]; 
cfg.boundmethod='ci';
cfg.cilevel = .95;
ft_boundedline_nightly(cfg,[],freq_norm_pos_test_ga,freq_norm_neg_test_ga)
fig1.Children.YLim(1)=0;
hold on
fill([foi(1) foi(1) foi(2) foi(2)],[0 max(fig1.Children.YLim) max(fig1.Children.YLim) 0],[216/255 216/255 216/255],'FaceAlpha',.5,'LineStyle','none')
legend({'Survivor','Non-survivor'})
title('All Chans')
% plot the topographies (plot 2/2)
subplot(143)
cfg=[];
cfg.marker = 'off';
cfg.comment = 'no';
cfg.xlim=foi;
cfg.zlim = [4.5e-3 9.5e-3];
cfg.layout='standard_1020.elc';
ft_topoplotER(cfg,freq_norm_pos_test_ga)
ylabel('Power (a.u.)')
xlabel('Frequency (Hz)')
colormap(flipud(brewermap(64,'RdBu')))
title(sprintf('Survivors\nfoi: %i - %i Hz',foi(1),foi(2)));
fig1_s3=subplot(144);
cfg=[];
cfg.marker = 'off';
cfg.comment = 'no';
cfg.xlim=foi;
cfg.zlim = [4.5e-3 9.5e-3];
cfg.layout='standard_1020.elc';
ft_topoplotER(cfg,freq_norm_neg_test_ga)
colormap(flipud(brewermap(64,'RdBu')))
title(sprintf('Non-survivors\nfoi: %i - %i Hz',foi(1),foi(2)));
fig1.Color = [1 1 1];
fig1.Position = [404 548 1203 430];
originalSize2 = fig1_s3.Position;
% hold on;
colorbar
fig1_s3.Position=originalSize2;
% export_fig('Z:\\USERS_DATA\\tkusterm\\resting_state\\Figures\\nothresh\\spectral\\prediction\\d1_spectra_10Hz_test_set_excluded_patients_%i.png','-painters','-r300')
%% 

%%%%%%%%%%% GROUP LEVEL COMPARISON POWER %%%%%%%%%%%

%% prepare neighbourhood structure 
cfg=[];
cfg.method         = 'template';%, 'triangulation' or 'template'
cfg.template       = 'gtec62_neighbours.mat';%name of the template fil
cfg.layout         = lay; %filename of the layout, see FT_PREPARE_LAYOUT
nghbrs             = ft_prepare_neighbours(cfg);
%% compute statistic for a comparison of outcome groups

% uncomment one of following depending on groups to compare
% training set
dat_a_stat = freq_all_norm_pos;
dat_b_stat = freq_all_norm_neg;
% test set
% dat_a_stat = freq_all_norm_pos_test;
% dat_b_stat = freq_all_norm_neg_test;

foi= [1 40];
cfg = [];
elecs = p_layout('ladybird');         
cfg.channel          = {'EEG'};
cfg.frequency        = foi;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
cfg.neighbours       = nghbrs;
cfg.design           = [ones(1,length(dat_a_stat)),ones(1,length(dat_b_stat))*2];
cfg.ivar             = 1;
[stat] = ft_freqstatistics(cfg, dat_a_stat{:}, dat_b_stat{:});
%% get grand average w/o retaining individuals
dat_a = ft_freqgrandaverage([],dat_a_stat{:});
dat_b = ft_freqgrandaverage([],dat_b_stat{:});
% keep individuals for boundedline plot
cfg = [];
cfg.keepindividual = 'yes';
dat_a_indiv_ga = ft_freqgrandaverage(cfg,dat_a_stat{:});
dat_b_indiv_ga = ft_freqgrandaverage(cfg,dat_b_stat{:});
% select the foi
cfg= [];
cfg.frequency = foi;
dat_a=ft_selectdata(cfg,dat_a);
dat_b=ft_selectdata(cfg,dat_b);
dat_a.mask=stat.mask;
dat_b.mask=stat.mask;
%% plot possible differences
% execute one of the combinations desired for plotting

% get time frame of cluster
clear idx_cluster idx_cluster_2
% idx_cluster(1)=min(stat.freq(any(stat.mask)));
% idx_cluster(2)=max(stat.freq(any(stat.mask)));

% % if no cluster exist, alternatively use n.s. ones to plot their extent
% idx_cluster(1)=min(stat.freq(any(stat.posclusterslabelmat==1)));
% idx_cluster(2)=max(stat.freq(any(stat.posclusterslabelmat==1)));
% idx_cluster(1)=min(stat.freq(any(stat.negclusterslabelmat==1)));
% idx_cluster(2)=max(stat.freq(any(stat.negclusterslabelmat==1)));


% multiple clusters (one pos and one neg here, e.g. for day 1 comp of surv and non-survs)
% idx_cluster(1)=min(stat.freq(any(stat.posclusterslabelmat==1)));
% idx_cluster(2)=max(stat.freq(any(stat.posclusterslabelmat==1)));
% idx_cluster_2(1)=min(stat.freq(any(stat.negclusterslabelmat==1)));
% idx_cluster_2(2)=max(stat.freq(any(stat.negclusterslabelmat==1)));
%% Plot spectra for paper
fig1=figure;
% set the cluster extent via grey bars at bottom of plot
fill([idx_cluster(1) idx_cluster(1) idx_cluster(2) idx_cluster(2)],[0 0.003 0.003 0],[166/255 166/255 166/255],'FaceAlpha',.5,'LineStyle','None');
if exist('idx_cluster_2','var')
    hold on;fill([idx_cluster_2(1) idx_cluster_2(1) idx_cluster_2(2) idx_cluster_2(2)],[0 0.003 0.003 0],[166/255 166/255 166/255],'FaceAlpha',.5,'LineStyle','None');
end
% plot the lines with confidence intervals
cfg = [];
cfg.channel = {stat.label{any(stat.posclusterslabelmat==1,2)}}; % plot significant channels
cfg.layout = lay;
cfg.ylim = [0 0.05];
cfg.boundmethod='ci';
cfg.cilevel = .95;
cfg2.frequency=[2 40];
ft_boundedline_nightly(cfg,[],dat_a_indiv_ga,dat_b_indiv_ga)
% change appearance
allines=findall(fig1,'Type','line');
set(allines,'LineWidth',2)
allaxes=findall(fig1,'Type','axes');
set(allaxes,'LineWidth',2)
hold on;
fig_tmp=gca;
hold on
% add legend
% [leg,icons]=legend({sprintf('Training ({\\itn}=%i)',numel(plot_a(~isnan(plot_a)))+numel(plot_b(~isnan(plot_b)))),sprintf('Test ({\\itn}=%i)',numel(plot_c(~isnan(plot_c)))+numel(plot_d(~isnan(plot_d))))},'Location','SouthEast','FontSize',11);
[leg,icon]=legend({'favorable outcome','unfavorable outcome'},'FontSize',12);
% leg.Position=[0.7246    0.1632    0.1254    0.0782];

% enlarge points in figure and legend
% set(findall(gca),'LineWidth',2);
% set(findobj(gca,'Type','Scatter'),'SizeData',400);
% set(findall(gca,'Type','Text'),'FontSize',11);
pause(1)
% icon(3).Children.MarkerSize=25;
% icon(4).Children.MarkerSize=25;
%
% [leg,icon]=legend({'Survivors','Non-survivors'});
xlabel('Frequency (Hz)')
ylabel('Normalized power')
% change color of lines
col1=[0/255 188/255 213/255];
col2=[255/255 94/255 105/255];
fig1.Children(2).Children(2).Color = col1;
fig1.Children(2).Children(4).FaceColor = col1;
fig1.Children(2).Children(1).Color = col2;
fig1.Children(2).Children(3).FaceColor = col2;
fig1.Position = [680 708 480 270];
fig1.Color=[1 1 1];
fig1.PaperPositionMode='auto';
% change aesthetics of legend
    % change color of legend markers
pause(1)
icon(3).FaceColor=col1;
icon(4).FaceColor=col2;
set(gca, 'box', 'off')
% change font sizes
allaxes=findall(fig1,'Type','axes');
set(allaxes,'FontSize',12)
allaxes=findall(fig1,'Type','legend');
set(allaxes,'FontSize',12)
%poscluster
text(9,.018,sprintf('{\\itp} = %.3f',stat.posclusters(1).prob),'FontSize',12)
%negcluster
text(32.5,.008,sprintf('{\\itp} = %.3f',stat.negclusters(1).prob),'FontSize',12)
fig1.Position=[834   390   885   388];
%%% PNG
% training
% print('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_spectra-training_1_40Hz_squareshape.png','-painters','-r600','-dpng')
% test
% print('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_spectra-test_inc_excluded_inc_comorbid_1_40Hz_squareshape.png','-painters','-r600','-dpng')
%%% PDF
% training
% print('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_spectra-training_1_40Hz.pdf','-painters','-dpdf')
% test
% print('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_spectra-test_inc_excluded_inc_comorbid_1_40Hz.pdf','-painters','-dpdf')
%% plot topographies for survivors, non-survivors and t-map of differences between the two in cluster
elecs=p_layout('ladybird'); 
fig1=figure;
cfg = [];
cfg.xlim = idx_cluster;
% cfg.zlim= [4.4E-3 8.7E-3];
cfg.zlim= [0.006 0.012]; % rom training set on day 1
cfg.marker='off';
cfg.comment='no';
cfg.highlight = 'on';
cfg.highlightchannel = {stat.label{any(stat.posclusterslabelmat,2)}};
cfg.highlightsymbol='.';
cfg.highlightsize = 6;
cfg.layout=lay;
% plot first group topograph
subplot(1,3,1);
ft_topoplotER(cfg,dat_a)
colormap(flipud(brewermap(64,'RdBu')))
% plot second group topograph
subplot(1,3,2);
ft_topoplotER(cfg,dat_b)
% plot t-statistic of difference in pos vs. neg outcome
subplot(1,3,3);
cfg.parameter='stat';
cfg.zlim=[0 3];
ft_topoplotER(cfg,stat)
oPos1=fig1.Children(1).Position;
fig1.Color=[1 1 1];
fig1.PaperPositionMode='auto';

% export_fig('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_topo-training.png','-r1000')
%% plot legends so they can be used/manipulated outside of MATLAB
fig1=figure;
% colorbar for normalized power
subplot(1,2,1)
colormap(flipud(brewermap(64,'RdBu')))
cb=colorbar('southoutside');
% modify colorbar
cb.Label.String = 'normalized power';
cb.Label.FontSize = 12;
cb.Ticks = [0.0050    0.0060    0.0070    0.0080];
cb.Ticks = [0.016    0.02    0.024];
caxis([0.016 0.024])
cb.Ticks = [0.006    0.008  0.010  0.012];
caxis([0.006 0.012])

% colorbar for statistical values (t-statistic)
subplot(1,2,2)
colormap(flipud(brewermap(64,'RdBu')))
cb=colorbar('southoutside');
% modify colorbar
cb.Label.String = 't-value';
cb.Label.FontSize = 12;
cb.Ticks = [0 1 2 3];
caxis([0 3])
% figure position
fig1.Position = [462 558 778 420];
fig1.Color = [1 1 1];
% set linewidth and fontsize of legend
allaxes=findall(fig1,'Type','colorbar');
set(allaxes,'LineWidth',1.5)
set(allaxes,'FontSize',12)
% PDF
% export_fig('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_legends-training.pdf')
% PNG
% export_fig('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_legends-training.png','-r1000')
%% combined plot of spectra, topographies and multiplot
%%% multiplot
fig1=figure;
subplot(2,4,[3 4 7 8]);
stat.posmask = stat.posclusterslabelmat==1;
stat.negmask = stat.negclusterslabelmat==1;
cfg = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskparameter = 'posmask';
cfg.layout = lay;
ft_multiplotER(cfg,stat)
fig1.Position=[119 351 1569 627];
%%% spectrum
cfg = [];
% cfg.parameter = 'powspctrm';
% elecs=p_layout('ladybird'); = {stat.label{any(stat.mask,2)}}; % plot significant channels
elecs=p_layout('ladybird'); 
cfg.channel = {stat.label{any(stat.posclusterslabelmat==1,2)}}; % plot significant channels
% cfg.maskparameter ='mask';
cfg.layout = lay;
subplot(2,4,[1 2]);
cfg.ylim = [0 0.05];
cfg.frequency=[1 40];
cfg.boundmethod='ci';
cfg.cilevel = .95;
ft_boundedline(cfg,[],dat_a_indiv_ga,dat_b_indiv_ga)
alllines=findall(fig1,'Type','line');
set(alllines,'LineWidth',1.5)
allaxes=findall(fig1,'Type','axes');
set(allaxes,'LineWidth',1.5)
hold on;
fig_tmp=gca;
hold on
% % plot the cluster extent on spectral lines either (in order of code): 
% % a) fill a rectangle with grey around frequencies
% % b) dashed line rectangle
% % c) add black bar at top of plot
% % d) shade area between curves (multiple lines)
fill([idx_cluster(1) idx_cluster(1) idx_cluster(2) idx_cluster(2)],[fig_tmp.YLim(1) fig_tmp.YLim(2) fig_tmp.YLim(2) fig_tmp.YLim(1)],[166/255 166/255 166/255],'FaceAlpha',.6,'LineStyle','None')
% rectangle('Position',[idx_cluster(1) fig_tmp.YLim(1) idx_cluster(2)-idx_cluster(1) fig_tmp.YLim(2)-fig_tmp.YLim(1)],'LineStyle','--')
% line([idx_cluster(1) idx_cluster(2)],[fig_tmp.YLim(2) fig_tmp.YLim(2)],'LineWidth',5,'Color',[0 0 0])
% % get indices and shade area between curves
% min_freq_idx=find(freq_norm_pos_ga.freq==idx_cluster(1));
% max_freq_idx=find(freq_norm_pos_ga.freq==idx_cluster(2));
% x2 = [freq_norm_pos_ga.freq(min_freq_idx:max_freq_idx), fliplr(freq_norm_pos_ga.freq(min_freq_idx:max_freq_idx))];
% inBetween = [squeeze(mean(mean(freq_norm_pos_ga.powspctrm(:,:,min_freq_idx:max_freq_idx),2),1))', fliplr(squeeze(mean(mean(freq_norm_neg_ga.powspctrm(:,:,min_freq_idx:max_freq_idx),2),1))')];
% fill(x2, inBetween, [166/255 166/255 166/255],'FaceAlpha',.6,'LineStyle','None');

%%% add legend and labels
legend({'Survivors','Non-survivors'})
xlabel('Frequency (Hz)')
ylabel('Normalized power')
% change line colors
col1=[0/255 188/255 213/255];
col2=[255/255 94/255 105/255];
fig1.Children(3).Children(2).Color = col1;
fig1.Children(3).Children(4).FaceColor = col1;
fig1.Children(3).Children(3).Color = col2;
fig1.Children(3).Children(5).FaceColor = col2;

%%% plot topography
cfg = [];
cfg.xlim = idx_cluster;
cfg.marker='off';
cfg.comment='no';
cfg.parameter = 'stat';
cfg.layout=lay;
plcfg.zlim=[0 3.2];
subplot(2,4,[5 6]);
ft_topoplotER(cfg,stat)
colormap(flipud(brewermap(64,'RdBu')))
title(sprintf('freq: %.2f - %.2f Hz',idx_cluster(1),idx_cluster(2)))
fig1.Color=[1 1 1];
fig1.PaperPositionMode='auto';
cb=colorbar;
cb.Label.String = 't-value';
cb.Label.FontSize = 12;
text(5,-.45,'t-values','FontSize',12)
% export_fig('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_stat_normalized-overlap.png','-painters','-r300')
%%


%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Predictions %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% select frequencies for which to predict outcome (e.g. from extent of cluster)
cfg = [];
% cfg.frequency = [10 10];
cfg.frequency = [5.2000 13.2000];
cfg.avgoverfreq = 'yes';
cfg.channel = 'EEG';
cfg.avgoverchan = 'yes';
% train
train_frq_pos = ft_selectdata(cfg,freq_norm_pos_ga);
train_frq_neg = ft_selectdata(cfg,freq_norm_neg_ga);
% test
test_frq_pos = ft_selectdata(cfg,freq_norm_pos_test_ga);
test_frq_neg = ft_selectdata(cfg,freq_norm_neg_test_ga);
%% find best ppv
adjust_prev='no'; % whether adjustment for prevalence in overall population should be applied
prev=0.5;% prevalence
ci_p_val=0.95;
z_val=norminv(ci_p_val);

% % 10 Hz power
dat_a=train_frq_pos.powspctrm;
dat_b=train_frq_neg.powspctrm;
dat_c=test_frq_pos.powspctrm;
dat_d=test_frq_neg.powspctrm;
title_1 = '10 Hz';
true_pos = 'above'; % whether subjects 'above' or 'below' the threshold should be considered to have predicted good outcome
%% compute ppv for all pruning levels 
TN=[];FN=[];FP=[];TP=[];
ppv=[];npv=[];sens=[];spec=[];acc=[];
ppv_ci=[];npv_ci=[];sens_ci=[];spec_ci=[];acc_ci=[];

thresh_pred=[];
% values for which to predict
thresh_pred=sort(unique([dat_a;dat_b]));
% compute prediction at each threshold value
for thresh_iter=1:numel(thresh_pred)
    if strcmpi(true_pos,'below')
        TN(thresh_iter)=length(find(dat_b' > thresh_pred(thresh_iter)));
        FN(thresh_iter)=length(find(dat_a' > thresh_pred(thresh_iter)));
        FP(thresh_iter)=length(find(dat_b' <= thresh_pred(thresh_iter)));
        TP(thresh_iter)=length(find(dat_a' <= thresh_pred(thresh_iter)));
    elseif strcmpi(true_pos,'above')
        TN(thresh_iter)=length(find(dat_b' < thresh_pred(thresh_iter)));
        FN(thresh_iter)=length(find(dat_a' < thresh_pred(thresh_iter)));
        FP(thresh_iter)=length(find(dat_b' >= thresh_pred(thresh_iter)));
        TP(thresh_iter)=length(find(dat_a' >= thresh_pred(thresh_iter)));
    end
    
    %
    TP_iter=TP(thresh_iter);
    FP_iter=FP(thresh_iter);
    TN_iter=TN(thresh_iter);
    FN_iter=FN(thresh_iter);
    
    % accuracy
    [acc(thresh_iter),acc_ci{thresh_iter}]=...
        binofit((TP_iter+TN_iter),(TP_iter+FP_iter+FN_iter+TN_iter),1-ci_p_val);
    
    % sensitivity
    [sens(thresh_iter),sens_ci{thresh_iter}]=...
        binofit(TP_iter,TP_iter+FN_iter,1-ci_p_val);
    
    % specificity
    [spec(thresh_iter),spec_ci{thresh_iter}]=...
        binofit(TN_iter,TN_iter+FP_iter,1-ci_p_val);
    
    % bifurcation for whether an adjustment for prevalence in population should be made for PPV and NPV
    if strcmp(adjust_prev,'no')
        % positive predictive value (PPV)
        [ppv(thresh_iter),ppv_ci{thresh_iter}]=...
            binofit(TP_iter,TP_iter+FP_iter,1-ci_p_val);
        % negative predictive value (NPV)
        [npv(thresh_iter),npv_ci{thresh_iter}]=...
            binofit(TN_iter,TN_iter+FN_iter);
    elseif strcmp(adjust_prev,'yes') % adjust for actual prevalence numbers
        % positive predictive value (PPV)
        ppv(thresh_iter)=(sens(thresh_iter)*prev)/...
            (sens(thresh_iter)*prev+(1-spec(thresh_iter))*(1-prev));
        ppv_tmp = ppv(thresh_iter);
        ppv_se=sqrt((ppv_tmp*(1-ppv_tmp))/(TP_iter+FP_iter));
        ppv_ci{thresh_iter}=[ppv_tmp-z_val*ppv_se, ppv_tmp+z_val*ppv_se];
        % negative predictive value (NPV)
        npv(thresh_iter)=(spec(thresh_iter)*(1-prev))/...
            (spec(thresh_iter)*(1-prev)+(1-sens(thresh_iter))*prev);
        npv_tmp=npv(thresh_iter);
        npv_se=sqrt((npv_tmp*(1-npv_tmp))/(TN_iter+FN_iter));
        npv_ci{thresh_iter}=[npv_tmp-z_val*npv_se, npv_tmp+z_val*npv_se];
    end
end
%% identify highest ppv across pruning levels and thresholds

% find best threshold per pruning level with min 10 predictions in PPV
pred_vals=[TP;FP;TN;FN;...
    ppv;npv;sens;spec;acc;thresh_pred']';
% find maximal ppv with min 10 correct TP
tmp_1=pred_vals;
tmp_1(TP<10,:)=0;
[~,idx]=max(tmp_1(:,5)); % finds max POSITIVE PREDICTIVE VALUE
% [~,idx]=max(tmp_1(:,9)); % finds max ACCURACY
max_ppv=tmp_1(idx,:);
thresh_set = max_ppv(10); % value to use for prediction in test set
% thresh_set = 0.0092; % original PPV tuned value from original training set
%% compute prediction values in the test set
% thresh_num=100;
prev=0.5;% prevalence
adjust_prev='no'; % adjustment for prevalence in population: 'yes' or 'no'
ci_p_val=0.95;
z_val=norminv(ci_p_val);

TN=[];FN=[];FP=[];TP=[];
ppv=[];npv=[];sens=[];spec=[];acc=[];
ppv_ci=[];npv_ci=[];sens_ci=[];spec_ci=[];acc_ci=[];
for thresh_iter=1:numel(thresh_set)    
    if strcmpi(true_pos,'below')
        % training set
        TN_train(thresh_iter)=length(find(dat_b' > thresh_set(thresh_iter)));
        FN_train(thresh_iter)=length(find(dat_a' > thresh_set(thresh_iter)));
        FP_train(thresh_iter)=length(find(dat_b' <= thresh_set(thresh_iter)));
        TP_train(thresh_iter)=length(find(dat_a' <= thresh_set(thresh_iter)));
        % test set        
        TN_test(thresh_iter)=length(find(dat_d' > thresh_set(thresh_iter)));
        FN_test(thresh_iter)=length(find(dat_c' > thresh_set(thresh_iter)));
        FP_test(thresh_iter)=length(find(dat_d' <= thresh_set(thresh_iter)));
        TP_test(thresh_iter)=length(find(dat_c' <= thresh_set(thresh_iter)));
    elseif strcmpi(true_pos,'above')
        % training set
        TN_train(thresh_iter)=length(find(dat_b' < thresh_set(thresh_iter)));
        FN_train(thresh_iter)=length(find(dat_a' < thresh_set(thresh_iter)));
        FP_train(thresh_iter)=length(find(dat_b' >= thresh_set(thresh_iter)));
        TP_train(thresh_iter)=length(find(dat_a' >= thresh_set(thresh_iter)));
        % test set        
        TN_test(thresh_iter)=length(find(dat_d' < thresh_set(thresh_iter)));
        FN_test(thresh_iter)=length(find(dat_c' < thresh_set(thresh_iter)));
        FP_test(thresh_iter)=length(find(dat_d' >= thresh_set(thresh_iter)));
        TP_test(thresh_iter)=length(find(dat_c' >= thresh_set(thresh_iter)));
    end
    % training
    [acc_train(thresh_iter),acc_ci_train{thresh_iter}]=binofit((TP_train(thresh_iter)+TN_train(thresh_iter)),(TP_train(thresh_iter)+FP_train(thresh_iter)+FN_train(thresh_iter)+TN_train(thresh_iter)),0.05);
    [sens_train(thresh_iter),sens_ci_train{thresh_iter}]=binofit(TP_train(thresh_iter),TP_train(thresh_iter)+FN_train(thresh_iter),0.05);
    [spec_train(thresh_iter),spec_ci_train{thresh_iter}]=binofit(TN_train(thresh_iter),TN_train(thresh_iter)+FP_train(thresh_iter),0.05);
    % test
    [acc_test(thresh_iter),acc_ci_test{thresh_iter}]=binofit((TP_test(thresh_iter)+TN_test(thresh_iter)),(TP_test(thresh_iter)+FP_test(thresh_iter)+FN_test(thresh_iter)+TN_test(thresh_iter)),0.05);
    [sens_test(thresh_iter),sens_ci_test{thresh_iter}]=binofit(TP_test(thresh_iter),TP_test(thresh_iter)+FN_test(thresh_iter),0.05);
    [spec_test(thresh_iter),spec_ci_test{thresh_iter}]=binofit(TN_test(thresh_iter),TN_test(thresh_iter)+FP_test(thresh_iter),0.05);
    
    % training set - bifurcation whether to adjust for the prevalence defined above
    if strcmp(adjust_prev,'no')
        [ppv_train(thresh_iter),ppv_ci_train{thresh_iter}]=binofit(TP_train(thresh_iter),TP_train(thresh_iter)+FP_train(thresh_iter),0.05);
        [npv_train(thresh_iter),npv_ci_train{thresh_iter}]=binofit(TN_train(thresh_iter),TN_train(thresh_iter)+FN_train(thresh_iter),0.05);
    elseif strcmp(adjust_prev,'yes') % adjust for actual prevalence numbers
        % positive predictive value (PPV)
        ppv_train(thresh_iter)=(sens_train(thresh_iter)*prev)/...
            (sens_train(thresh_iter)*prev+(1-spec_train(thresh_iter))*(1-prev));
        ppv_tmp_train = ppv_train(thresh_iter);
        ppv_se_train=sqrt((ppv_tmp_train*(1-ppv_tmp_train))/(TP_train(thresh_iter)+FP_train(thresh_iter)));
        ppv_ci_train{thresh_iter}=[ppv_tmp_train-z_val*ppv_se_train, ppv_tmp_train+z_val*ppv_se_train];
        % negative predictive value (NPV)
        npv_train(thresh_iter)=(spec_train(thresh_iter)*(1-prev))/...
            ((1-sens_train(thresh_iter))*prev+spec_train(thresh_iter)*(1-prev));
        npv_tmp_train=npv_train(thresh_iter);
        npv_se_train=sqrt((npv_tmp_train*(1-npv_tmp_train))/(TN_train(thresh_iter)+FN_train(thresh_iter)));
        npv_ci_train{thresh_iter}=[npv_tmp_train-z_val*npv_se_train, npv_tmp_train+z_val*npv_se_train];
    end
    % test - bifurcation whether to adjust for the prevalence defined above
    if strcmp(adjust_prev,'no')
        [ppv_test(thresh_iter),ppv_ci_test{thresh_iter}]=binofit(TP_test(thresh_iter),TP_test(thresh_iter)+FP_test(thresh_iter),0.05);
        [npv_test(thresh_iter),npv_ci_test{thresh_iter}]=binofit(TN_test(thresh_iter),TN_test(thresh_iter)+FN_test(thresh_iter),0.05);
    elseif strcmp(adjust_prev,'yes') % adjust for actual prevalence numbers
        % positive predictive value (PPV)
        ppv_test(thresh_iter)=(sens_test(thresh_iter)*prev)/...
            (sens_test(thresh_iter)*prev+(1-spec_test(thresh_iter))*(1-prev));
        ppv_tmp_test = ppv_test(thresh_iter);
        ppv_se_test=sqrt((ppv_tmp_test*(1-ppv_tmp_test))/(TP_test(thresh_iter)+FP_test(thresh_iter)));
        ppv_ci_test{thresh_iter}=[ppv_tmp_test-z_val*ppv_se_test, ppv_tmp_test+z_val*ppv_se_test];
        % negative predictive value (NPV)
        npv_test(thresh_iter)=(spec_test(thresh_iter)*(1-prev))/...
            ((1-sens_test(thresh_iter))*prev+spec_test(thresh_iter)*(1-prev));
        npv_tmp_test=npv_test(thresh_iter);
        npv_se_test=sqrt((npv_tmp_test*(1-npv_tmp_test))/(TN_test(thresh_iter)+FN_test(thresh_iter)));
        npv_ci_test{thresh_iter}=[npv_tmp_test-z_val*npv_se_test, npv_tmp_test+z_val*npv_se_test];
    end;    
end;
%%
% move to plotting variables
plot_a=dat_a';
plot_b=dat_b';
plot_c=dat_c';
plot_d=dat_d';

% add jitter in plotting
rand_a=rand(size(plot_a));
rand_b=rand(size(plot_b));
rand_c=rand(size(plot_c));
rand_d=rand(size(plot_d));
%% PLOT - adaptive threshold plotting

alpha_deg = 1; % opacity
random_degree=0.2; % magnitude of horizontal scatter jitter

% the color values were (probably) derived from the color scheme of gramm (itself using the ggplot scheme I think) and creating a more
% "red" version plus a lighter red version for differentiation - same for blue colors
colors_scatter={[0 .57 .84],[1 .37 .41],[0 .74 .84],[255/255 102/255 204/255]};% dark blue, dark red,light blue, purple/red;

idx=1;

grps = {'favorable','unfavorable'};
% title_plot = sprintf('%s - %.1f%% pruned',title_1,100-thresh_to_sel);
% prepare aesthetics (e.g. boxes around individual points)
min_box = min([plot_a,plot_b,plot_c,plot_d])-range([plot_a,plot_b,plot_c,plot_d])*0.05;
thresh_low = thresh_set - min_box - range([plot_a,plot_b,plot_c,plot_d])*0.005;
thresh_high = thresh_set + range([plot_a,plot_b,plot_c,plot_d])*0.005;
max_box = max([plot_a,plot_b,plot_c,plot_d])+range([plot_a,plot_b,plot_c,plot_d])*0.05 - thresh_high;
pred_text=max_box+thresh_high-range([plot_a,plot_b,plot_c,plot_d])*(0.05:0.07:0.40);
pred_text_2 = pred_text(end)-range([plot_a,plot_b,plot_c,plot_d])*(0.10:0.07:0.45);
thresh_text = thresh_set + range([plot_a,plot_b,plot_c,plot_d])*0.05;

label_low = min_box+thresh_low/2;
label_high = thresh_high+max_box/2;

fig1=figure;hold on;
% plot survivors/non-survivors by training/test set
scatter(-0.05+1-rand_a*random_degree,plot_a,'Marker','.','MarkerEdgeColor',[0 .57 .84],'MarkerFaceColor',[0 .57 .84],'MarkerEdgeAlpha',alpha_deg)
scatter(0.05+1+rand_c*random_degree,plot_c,'Marker','.','MarkerEdgeColor',[1 .37 .41],'MarkerFaceColor',[1 .37 .41],'MarkerEdgeAlpha',alpha_deg)
scatter(-0.05+2-rand_b*random_degree,plot_b,'Marker','.','MarkerEdgeColor',[0 .57 .84],'MarkerFaceColor',[0 .57 .84],'MarkerEdgeAlpha',alpha_deg)
scatter(0.05+2+rand_d*random_degree,plot_d,'Marker','.','MarkerEdgeColor',[1 .37 .41],'MarkerFaceColor',[1 .37 .41],'MarkerEdgeAlpha',alpha_deg)


ylabel('normalized power')
xlim([-.5 4.5]) % set x axis to create space to fit annotations
plot([0.3 2.7],[thresh_set thresh_set],'LineStyle','--','Color',[169/255 169/255 169/255]) % threshold line
% label columns (x-axis)
fig1.Children.XTick = [1 2];
fig1.Children.XTickLabel = grps;
h=xlabel('Outcome');
h.Position= [1.5 -0.0011 -1];
% rectangles around subjects
rectangle('Position', [0.7 min_box 0.6 thresh_low],'LineStyle','--','LineWidth',1.25)
rectangle('Position', [0.7 thresh_high 0.6 max_box],'LineStyle','--','LineWidth',1.25)
rectangle('Position', [1.7 min_box 0.6 thresh_low],'LineStyle','--','LineWidth',1.25)
rectangle('Position', [1.7 thresh_high 0.6 max_box],'LineStyle','--','LineWidth',1.25)

% gather prediction subject counts
if strcmp(true_pos,'above')
    low_s = [FN_train,FN_test];
    high_s = [TP_train,TP_test];
    high_ns = [FP_train,FP_test];
    low_ns = [TN_train,TN_test];
elseif strcmp(true_pos,'below')
    low_s = [TP_train,TP_test];
    high_s = [FN_train,FN_test];
    high_ns = [TN_train,TN_test];
    low_ns = [FP_train,FP_test];
end
text(-.15,label_high,sprintf('Training/Test \nn = (%s/%s)',num2str(high_s(1)),num2str(high_s(2))),'FontSize',11)
text(.1,label_low,['n = (',num2str(low_s(1)),'/',num2str(low_s(2)),')'],'FontSize',11)
text(2.35,label_low,['n = (',num2str(low_ns(1)),'/',num2str(low_ns(2)),')'],'FontSize',11)
text(2.35,label_high,['n = (',num2str(high_ns(1)),'/',num2str(high_ns(2)),')'],'FontSize',11)

% add prediction values for training set
text(3,pred_text(1),sprintf('{\\bfTraining} ({\\itn}=%i)',numel(plot_a)+numel(plot_b)))
text(3,pred_text(2),sprintf('PPV: %.2f (CI: %.2f-%.2f)',ppv_train,ppv_ci_train{idx}(1),ppv_ci_train{idx}(2)))
text(3,pred_text(3),sprintf('NPV: %.2f (CI: %.2f-%.2f)',npv_train,npv_ci_train{idx}(1),npv_ci_train{idx}(2)))
text(3,pred_text(4),sprintf('Sensitivity: %.2f (CI: %.2f-%.2f)',sens_train(idx),sens_ci_train{idx}(1),sens_ci_train{idx}(2)))
text(3,pred_text(5),sprintf('Specificity: %.2f (CI: %.2f-%.2f)',spec_train(idx),spec_ci_train{idx}(1),spec_ci_train{idx}(2)))
text(3,pred_text(6),sprintf('Accuracy: %.2f (CI: %.2f-%.2f)',acc_train(idx),acc_ci_train{idx}(1),acc_ci_train{idx}(2)))
text(-.3,thresh_text,sprintf('Threshold: %.5f',thresh_set))
% add prediction values for test set
text(3,pred_text_2(1),sprintf('{\\bfTest} ({\\itn}=%i)',numel(plot_c)+numel(plot_d)))
text(3,pred_text_2(2),sprintf('PPV: %.2f (CI: %.2f-%.2f)',ppv_test(idx),ppv_ci_test{idx}(1),ppv_ci_test{idx}(2)))
text(3,pred_text_2(3),sprintf('NPV: %.2f (CI: %.2f-%.2f)',npv_test(idx),npv_ci_test{idx}(1),npv_ci_test{idx}(2)))
text(3,pred_text_2(4),sprintf('Sensitivity: %.2f (CI: %.2f-%.2f)',sens_test(idx),sens_ci_test{idx}(1),sens_ci_test{idx}(2)))
text(3,pred_text_2(5),sprintf('Specificity: %.2f (CI: %.2f-%.2f)',spec_test(idx),spec_ci_test{idx}(1),spec_ci_test{idx}(2)))
text(3,pred_text_2(6),sprintf('Accuracy: %.2f (CI: %.2f-%.2f)',acc_test(idx),acc_ci_test{idx}(1),acc_ci_test{idx}(2)))

% change figure aesthetics
fig1.Color = [1 1 1];
fig1.Position=[834   390   885   588];
fig1.Children.FontSize=12;
% title(title_plot);

% add legend
[leg,icons]=legend({sprintf('Training ({\\itn}=%i)',numel(plot_a(~isnan(plot_a)))+numel(plot_b(~isnan(plot_b)))),
    sprintf('Test ({\\itn}=%i)',numel(plot_c(~isnan(plot_c)))+numel(plot_d(~isnan(plot_d))))},'Location','SouthEast','FontSize',12);
leg.Position=[0.7246    0.1632    0.1254    0.0782];

% enlarge points in figure and legend
set(findall(gca),'LineWidth',2);
set(findobj(gca,'Type','Scatter'),'SizeData',400);
set(findall(gca,'Type','Text'),'FontSize',11);
pause(1)
icons(3).Children.MarkerSize=25;
icons(4).Children.MarkerSize=25;

fig1.Children(2).YTickLabel(fig1.Children(2).YTick < 0)={''}; % remove ticks below zero
fig1.PaperPositionMode = 'auto';
%%
% PNG
% export_fig('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_prediction_spectra_5.2-13.2Hz_train_test_incl_excluded.png','-painters','-r300')
% PDF
% export_fig('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_prediction_spectra_5.2-13.2Hz_train_test_incl_excluded_incl_comorbid.pdf')
%%
% print('D:\\Google Drive\\Arbeit\\Lausanne\\resting_state\\Paper_Spectra\\Figures\\Spectra\\d1_prediction_spectra_5.2-13.2Hz_train_test_incl_excluded_inc_comorbid.pdf','-painters','-dpdf')
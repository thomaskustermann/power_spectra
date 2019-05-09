% script to export the predictions from power spectra and graph metrics to a text file
base_path = 'Z:\USERS_DATA\tkusterm\resting_state\DATA\RestingEEG\resting_data\';
%%
% day 1
training_pos=training_pos_d1;
training_neg=training_neg_d1;
test_pos=test_pos_d1;
test_neg=test_neg_d1;

% day 2
% training_pos=training_pos_d2;
% training_neg=training_neg_d2;
% test_pos=test_pos_d2;
% test_neg=test_neg_d2;

%% load PRE-saved graph metrics

file_to_load = 'graphm_2018_07_31_var_over_time.mat'; % thresholded - binarized - change over time (var) - failed thresholdings, fixed degree calc
graphm=[];
thresh_val = 50:-2.5:10;
% training - survivors
for subj_iter=1:numel(training_pos)
    load_path = [base_path,training_pos{subj_iter},'5s_segm_nothresh\connectivity\',file_to_load]; 
    display(['Loading from: ',load_path])
    % load the individual graph results
    load(load_path)
    % %    clustering coefficient
    graphm.cl_coef.avg_pos_auc(subj_iter) = graphm_indiv.cl_coef.avg_pos_auc;
    % %     path length/distance
    graphm.dist.avg_pos_auc(subj_iter) = graphm_indiv.dist.avg_pos_auc;
    graphm.dist.std_pos_auc(subj_iter) = graphm_indiv.dist.std_pos_auc;
    % %     participation coefficient and modularity
    graphm.part_coef.mod_pos_auc(subj_iter) = graphm_indiv.part_coef.mod_pos_auc;
    graphm.part_coef.part_pos_mean_auc(subj_iter) = graphm_indiv.part_coef.part_pos_mean_auc;
    graphm.part_coef.part_pos_std_auc(subj_iter) = graphm_indiv.part_coef.part_pos_std_auc;
    for thresh_iter=1:numel(thresh_val)
        % %         clustering coefficient
        graphm.cl_coef.pos{thresh_iter}{subj_iter} = graphm_indiv.cl_coef.pos{thresh_iter};
        graphm.cl_coef.avg_pos{thresh_iter}(subj_iter) = graphm_indiv.cl_coef.avg_pos(thresh_iter);
        % %         path length/distance
        graphm.dist.pos{thresh_iter}{subj_iter} = graphm_indiv.dist.pos{thresh_iter};
        graphm.dist.avg_pos{thresh_iter}(subj_iter) = graphm_indiv.dist.avg_pos{thresh_iter};
        graphm.dist.std_pos{thresh_iter}(subj_iter) = graphm_indiv.dist.std_pos{thresh_iter};
        % %         participation coefficient and modularity
        graphm.part_coef.louv_pos{thresh_iter}{subj_iter} = graphm_indiv.part_coef.louv_pos{thresh_iter};
        graphm.part_coef.mod_pos{thresh_iter}{subj_iter} = graphm_indiv.part_coef.mod_pos{thresh_iter};
        graphm.part_coef.part_pos{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos{thresh_iter};
        graphm.part_coef.mod_pos_avg{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg{thresh_iter};
        graphm.part_coef.mod_pos_avg_var{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg_var{thresh_iter};
        graphm.part_coef.part_pos_mean_iter{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos_mean_iter{thresh_iter};
        graphm.part_coef.part_pos_mean{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_mean{thresh_iter};
        graphm.part_coef.part_pos_std{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_std{thresh_iter};
        % %         node degree
        graphm.degree.pos{thresh_iter}{subj_iter} = graphm_indiv.degree.pos{thresh_iter};
        graphm.degree.avg_pos{thresh_iter}(subj_iter) = graphm_indiv.degree.avg_pos{thresh_iter};
        % %         threshold failures
        graphm.thresh_fail.pos{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos{thresh_iter};
        graphm.thresh_fail.pos_perc{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos_perc{thresh_iter};
    end
end

% training - non-survivors
for subj_iter=1:numel(training_neg)
    load_path = [base_path,training_neg{subj_iter},'5s_segm_nothresh\connectivity\',file_to_load]; 
    display(['Loading from: ',load_path])
    % load the individual graph results
    load(load_path)
    % %    clustering coefficient
    graphm.cl_coef.avg_neg_auc(subj_iter) = graphm_indiv.cl_coef.avg_pos_auc;
    % %     path length/distance
    graphm.dist.avg_neg_auc(subj_iter) = graphm_indiv.dist.avg_pos_auc;
    graphm.dist.std_neg_auc(subj_iter) = graphm_indiv.dist.std_pos_auc;
    % %     participation coefficient and modularity
    graphm.part_coef.mod_neg_auc(subj_iter) = graphm_indiv.part_coef.mod_pos_auc;
    graphm.part_coef.part_neg_mean_auc(subj_iter) = graphm_indiv.part_coef.part_pos_mean_auc;
    graphm.part_coef.part_neg_std_auc(subj_iter) = graphm_indiv.part_coef.part_pos_std_auc;
    for thresh_iter=1:numel(thresh_val)
        % %         clustering coefficient
        graphm.cl_coef.neg{thresh_iter}{subj_iter} = graphm_indiv.cl_coef.pos{thresh_iter};
        graphm.cl_coef.avg_neg{thresh_iter}(subj_iter) = graphm_indiv.cl_coef.avg_pos(thresh_iter);
        % %         path length/distance
        graphm.dist.neg{thresh_iter}{subj_iter} = graphm_indiv.dist.pos{thresh_iter};
        graphm.dist.avg_neg{thresh_iter}(subj_iter) = graphm_indiv.dist.avg_pos{thresh_iter};
        graphm.dist.std_neg{thresh_iter}(subj_iter) = graphm_indiv.dist.std_pos{thresh_iter};
        % %         participation coefficient and modularity
        graphm.part_coef.louv_neg{thresh_iter}{subj_iter} = graphm_indiv.part_coef.louv_pos{thresh_iter};
        graphm.part_coef.mod_neg{thresh_iter}{subj_iter} = graphm_indiv.part_coef.mod_pos{thresh_iter};
        graphm.part_coef.part_neg{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos{thresh_iter};
        graphm.part_coef.mod_neg_avg{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg{thresh_iter};        
        graphm.part_coef.mod_neg_avg_var{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg_var{thresh_iter};
        graphm.part_coef.part_neg_mean_iter{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos_mean_iter{thresh_iter};
        graphm.part_coef.part_neg_mean{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_mean{thresh_iter};
        graphm.part_coef.part_neg_std{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_std{thresh_iter};
        % %         node degree
        graphm.degree.neg{thresh_iter}{subj_iter} = graphm_indiv.degree.pos{thresh_iter};
        graphm.degree.avg_neg{thresh_iter}(subj_iter) = graphm_indiv.degree.avg_pos{thresh_iter};
        % %         threshold failures
        graphm.thresh_fail.neg{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos{thresh_iter};
        graphm.thresh_fail.neg_perc{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos_perc{thresh_iter};
    end
end

% test - survivors
for subj_iter=1:numel(test_pos)
    load_path = [base_path,test_pos{subj_iter},'5s_segm_nothresh\connectivity\',file_to_load];
    display(['Loading from: ',load_path])
    % load the individual graph results
    load(load_path)
    % %    clustering coefficient
    graphm.cl_coef.avg_pos_auc_test(subj_iter) = graphm_indiv.cl_coef.avg_pos_auc;
    % %    path length/distance
    graphm.dist.avg_pos_auc_test(subj_iter) = graphm_indiv.dist.avg_pos_auc;
    graphm.dist.std_pos_auc_test(subj_iter) = graphm_indiv.dist.std_pos_auc;
    % %     participation coefficient and modularity
    graphm.part_coef.mod_pos_auc_test(subj_iter) = graphm_indiv.part_coef.mod_pos_auc;
    graphm.part_coef.part_pos_mean_auc_test(subj_iter) = graphm_indiv.part_coef.part_pos_mean_auc;
    graphm.part_coef.part_pos_std_auc_test(subj_iter) = graphm_indiv.part_coef.part_pos_std_auc;
    for thresh_iter=1:numel(thresh_val)
        % %         clustering coefficient
        graphm.cl_coef.pos_test{thresh_iter}{subj_iter} = graphm_indiv.cl_coef.pos{thresh_iter};
        graphm.cl_coef.avg_pos_test{thresh_iter}(subj_iter) = graphm_indiv.cl_coef.avg_pos(thresh_iter);
        % %         path length/distance
        graphm.dist.pos_test{thresh_iter}{subj_iter} = graphm_indiv.dist.pos{thresh_iter};
        graphm.dist.avg_pos_test{thresh_iter}(subj_iter) = graphm_indiv.dist.avg_pos{thresh_iter};
        graphm.dist.std_pos_test{thresh_iter}(subj_iter) = graphm_indiv.dist.std_pos{thresh_iter};
        % %         participation coefficient and modularity
        graphm.part_coef.louv_pos_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.louv_pos{thresh_iter};
        graphm.part_coef.mod_pos_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.mod_pos{thresh_iter};
        graphm.part_coef.part_pos_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos{thresh_iter};
        graphm.part_coef.mod_pos_avg_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg{thresh_iter};
        graphm.part_coef.mod_pos_avg_var_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg_var{thresh_iter};
        graphm.part_coef.part_pos_mean_iter_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos_mean_iter{thresh_iter};
        graphm.part_coef.part_pos_mean_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_mean{thresh_iter};
        graphm.part_coef.part_pos_std_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_std{thresh_iter};
        % %         node degree
        graphm.degree.pos_test{thresh_iter}{subj_iter} = graphm_indiv.degree.pos{thresh_iter};
        graphm.degree.avg_pos_test{thresh_iter}(subj_iter) = graphm_indiv.degree.avg_pos{thresh_iter};
        % %         threshold failures
        graphm.thresh_fail.pos_test{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos{thresh_iter};
        graphm.thresh_fail.pos_perc_test{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos_perc{thresh_iter};
    end
end

% test - non-survivors
for subj_iter=1:numel(test_neg)
    load_path = [base_path,test_neg{subj_iter},'5s_segm_nothresh\connectivity\',file_to_load];
    display(['Loading from: ',load_path])
    % load the individual graph results
    load(load_path)
    % %     clustering coefficient
    graphm.cl_coef.avg_neg_auc_test(subj_iter) = graphm_indiv.cl_coef.avg_pos_auc;
    % %     path length/distance
    graphm.dist.avg_neg_auc_test(subj_iter) = graphm_indiv.dist.avg_pos_auc;
    graphm.dist.std_neg_auc_test(subj_iter) = graphm_indiv.dist.std_pos_auc;
    % %     participation coefficient and modularity
    graphm.part_coef.mod_neg_auc_test(subj_iter) = graphm_indiv.part_coef.mod_pos_auc;
    graphm.part_coef.part_neg_mean_auc_test(subj_iter) = graphm_indiv.part_coef.part_pos_mean_auc;
    graphm.part_coef.part_neg_std_auc_test(subj_iter) = graphm_indiv.part_coef.part_pos_std_auc;
    for thresh_iter=1:numel(thresh_val)
        % %         clustering coefficient
        graphm.cl_coef.neg_test{thresh_iter}{subj_iter} = graphm_indiv.cl_coef.pos{thresh_iter};
        graphm.cl_coef.avg_neg_test{thresh_iter}(subj_iter) = graphm_indiv.cl_coef.avg_pos(thresh_iter);
        % %         path length/distance
        graphm.dist.neg_test{thresh_iter}{subj_iter} = graphm_indiv.dist.pos{thresh_iter};
        graphm.dist.avg_neg_test{thresh_iter}(subj_iter) = graphm_indiv.dist.avg_pos{thresh_iter};
        graphm.dist.std_neg_test{thresh_iter}(subj_iter) = graphm_indiv.dist.std_pos{thresh_iter};
        % %         participation coefficient and modularity
        graphm.part_coef.louv_neg_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.louv_pos{thresh_iter};
        graphm.part_coef.mod_neg_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.mod_pos{thresh_iter};
        graphm.part_coef.part_neg_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos{thresh_iter};
        graphm.part_coef.mod_neg_avg_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg{thresh_iter};
        graphm.part_coef.mod_neg_avg_var_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.mod_pos_avg_var{thresh_iter};
        graphm.part_coef.part_neg_mean_iter_test{thresh_iter}{subj_iter} = graphm_indiv.part_coef.part_pos_mean_iter{thresh_iter};
        graphm.part_coef.part_neg_mean_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_mean{thresh_iter};
        graphm.part_coef.part_neg_std_test{thresh_iter}(subj_iter) = graphm_indiv.part_coef.part_pos_std{thresh_iter};
        % %         node degree
        graphm.degree.neg_test{thresh_iter}{subj_iter} = graphm_indiv.degree.pos{thresh_iter};
        graphm.degree.avg_neg_test{thresh_iter}(subj_iter) = graphm_indiv.degree.avg_pos{thresh_iter};
        % %         threshhold failures
        graphm.thresh_fail.neg_test{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos{thresh_iter};
        graphm.thresh_fail.neg_perc_test{thresh_iter}(subj_iter) = graphm_indiv.thresh_fail.pos_perc{thresh_iter};
    end
end
%% add the frequency threshold information
thresh_freq = 0.01054; % value taken from 5.2 - 13.2 Hz day 1training set maximum PPV
%% load freq spectrum by outcome group
file_name = 'data_freq_mtm.mat';
% TRAINING
% pos outcome
for subj_iter = 1:numel(training_pos)
    file_to_load=[base_path,training_pos{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird'); 
    cfg.channel= {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_pos{subj_iter} = frq;
    % normalize for each subject
    freq_all_norm_pos{subj_iter} = freq_all_pos{subj_iter};
    freq_all_norm_pos{subj_iter}.powspctrm = bsxfun(@rdivide, freq_all_pos{subj_iter}.powspctrm, sum(freq_all_pos{subj_iter}.powspctrm,2));
end;
% neg outcome
for subj_iter = 1:numel(training_neg)
    file_to_load=[base_path,training_neg{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird'); 
    cfg.channel = {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_neg{subj_iter} = frq;
    % normalize for each subject
    freq_all_norm_neg{subj_iter} = freq_all_neg{subj_iter};
    freq_all_norm_neg{subj_iter}.powspctrm = bsxfun(@rdivide, freq_all_neg{subj_iter}.powspctrm, sum(freq_all_neg{subj_iter}.powspctrm,2));
end;

% TEST
% pos outcome
for subj_iter = 1:numel(test_pos)
    file_to_load=[base_path,test_pos{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird'); 
    cfg.channel= {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_pos_test{subj_iter} = frq;
    % normalize for each subject
    freq_all_norm_pos_test{subj_iter} = freq_all_pos_test{subj_iter};
    freq_all_norm_pos_test{subj_iter}.powspctrm = bsxfun(@rdivide, freq_all_pos_test{subj_iter}.powspctrm, sum(freq_all_pos_test{subj_iter}.powspctrm,2));
end;
% neg outcome
for subj_iter = 1:numel(test_neg)
    file_to_load=[base_path,test_neg{subj_iter},'\5s_segm_nothresh\freq\',file_name];
    load(file_to_load)
    display(file_to_load);
    cfg = [];
    cfg.frequency = [2 40];
    elecs=p_layout('ladybird'); 
    cfg.channel = {frq.label{~ismember(frq.label,{'EOGH','EOGV','ECG'})}}; % remove ECG/EOG channels from data
    frq=ft_selectdata(cfg,frq);
    freq_all_neg_test{subj_iter} = frq;
    % normalize for each subject
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
%% select and avg power spectra
cfg = [];
cfg.avgoverfreq = 'yes';
cfg.frequency = [5.2 13.2];
cfg.channel = 'EEG';
cfg.avgoverchan = 'yes';
% train
train_10_pos = ft_selectdata(cfg,freq_norm_pos_ga);
train_10_neg = ft_selectdata(cfg,freq_norm_neg_ga);
% test
test_10_pos = ft_selectdata(cfg,freq_norm_pos_test_ga);
test_10_neg = ft_selectdata(cfg,freq_norm_neg_test_ga);
%% check validity of power spectra predictions;
TP_train_frq = train_10_pos.powspctrm >= thresh_freq;
FP_train_frq = train_10_neg.powspctrm >= thresh_freq;
FN_train_frq = train_10_pos.powspctrm < thresh_freq;
TN_train_frq = train_10_neg.powspctrm < thresh_freq;

TP_test_frq = test_10_pos.powspctrm >= thresh_freq;
FP_test_frq = test_10_neg.powspctrm >= thresh_freq;
FN_test_frq = test_10_pos.powspctrm < thresh_freq;
TN_test_frq = test_10_neg.powspctrm < thresh_freq;

% add
% % code below used for value export
subj_a=training_pos;
subj_b=training_neg;
subj_c=test_pos;
subj_d=test_neg;

%
if size(test_pos,2) == 1
   test_pos = test_pos';
   test_neg = test_neg';
end

% variable gathering information for later export as text
thresh_all={};
% thresh_all(:,1)=['id';training_pos';training_neg';test_pos';test_neg'];
thresh_all(:,1)=['id';training_pos;training_neg;test_pos';test_neg'];
thresh_all(:,2)= ['grp_pred';...
    repmat({'training_surv'},numel(subj_a),1);...
    repmat({'training_nonsurv'},numel(subj_b),1);...
    repmat({'test_surv'},numel(subj_c),1);...
    repmat({'test_nonsurv'},numel(subj_d),1)];

title_1 = '10hz_power';
title1_low =lower(title_1);
dat_outcome=[train_10_pos.powspctrm;train_10_neg.powspctrm;test_10_pos.powspctrm;test_10_neg.powspctrm]>=thresh_freq;
dat_outcome = double(dat_outcome);
dat_outcome(isnan([train_10_pos.powspctrm;train_10_neg.powspctrm;test_10_pos.powspctrm;test_10_neg.powspctrm]))=NaN;
dat_val = [train_10_pos.powspctrm;train_10_neg.powspctrm;test_10_pos.powspctrm;test_10_neg.powspctrm];
dat_val(isnan([train_10_pos.powspctrm;train_10_neg.powspctrm;test_10_pos.powspctrm;test_10_neg.powspctrm]))=NaN;
thresh_all(:,7)=[title1_low;num2cell(dat_outcome)];
thresh_all(:,12)=[[title1_low,'_val'];num2cell(dat_val)];
%% Add the graph metrics
% select whether data should be log transformed to reign in some outliers
log_transform = 0;
if log_transform == 1
    warning('Log-transform not applied to values in "thresh_all"')
end

% % % clustering coefficient
thresh_to_sel=12.5;thresh_pred=0.007299;true_pos = 'above';
thresh_sel=find(thresh_val == thresh_to_sel); % pruning level to select
dat_a=graphm.cl_coef.avg_pos{thresh_sel};
dat_b=graphm.cl_coef.avg_neg{thresh_sel};
dat_c=graphm.cl_coef.avg_pos_test{thresh_sel};
dat_d=graphm.cl_coef.avg_neg_test{thresh_sel};
title_1 = 'Clustering coefficient';
% % code below used for value export
if strcmp(true_pos,'above')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']>=thresh_pred;
else strcmp(true_pos,'below')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']<=thresh_pred;
end
title1_low =lower(title_1);
dat_outcome = double(dat_outcome);
dat_outcome(isnan([dat_a';dat_b';dat_c';dat_d']))=NaN;
thresh_all(:,3)=[title1_low;num2cell(dat_outcome)];
thresh_all(:,8)=[[title1_low,' val'];num2cell([dat_a';dat_b';dat_c';dat_d'])];
thresh_all(:,13)=[[title1_low,'_fail_thresh_perc'];num2cell([graphm.thresh_fail.pos_perc{thresh_sel},graphm.thresh_fail.neg_perc{thresh_sel},...
    graphm.thresh_fail.pos_perc_test{thresh_sel},graphm.thresh_fail.neg_perc_test{thresh_sel}]')];



% % % modularity
thresh_to_sel=10;thresh_pred=0.000540;true_pos = 'above';
thresh_sel=find(thresh_val == thresh_to_sel); % pruning level to select
dat_a=graphm.part_coef.mod_pos_avg_var{thresh_sel};
dat_b=graphm.part_coef.mod_neg_avg_var{thresh_sel};
dat_c=graphm.part_coef.mod_pos_avg_var_test{thresh_sel};
dat_d=graphm.part_coef.mod_neg_avg_var_test{thresh_sel};
title_1 = 'Modularity';

% % code below used for value export
if strcmp(true_pos,'above')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']>=thresh_pred;
else strcmp(true_pos,'below')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']<=thresh_pred;
end
title1_low =lower(title_1);
dat_outcome = double(dat_outcome);
dat_outcome(isnan([dat_a';dat_b';dat_c';dat_d']))=NaN;
thresh_all(:,4)=[title1_low;num2cell(dat_outcome)];
thresh_all(:,9)=[[title1_low,' val'];num2cell([dat_a';dat_b';dat_c';dat_d'])];
thresh_all(:,14)=[[title1_low,'_fail_thresh_perc'];num2cell([graphm.thresh_fail.pos_perc{thresh_sel},graphm.thresh_fail.neg_perc{thresh_sel},...
    graphm.thresh_fail.pos_perc_test{thresh_sel},graphm.thresh_fail.neg_perc_test{thresh_sel}]')];


% % % path length/distance
thresh_to_sel=32.5;thresh_pred=0.000645;true_pos = 'above';
thresh_sel=find(thresh_val == thresh_to_sel); % pruning level to select
dat_a=graphm.dist.avg_pos{thresh_sel};
dat_b=graphm.dist.avg_neg{thresh_sel};
dat_c=graphm.dist.avg_pos_test{thresh_sel};
dat_d=graphm.dist.avg_neg_test{thresh_sel};
title_1 = 'Path length';

% code below used for value export
if strcmp(true_pos,'above')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']>=thresh_pred;
else strcmp(true_pos,'below')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']<=thresh_pred;
end
title1_low =lower(title_1);
dat_outcome = double(dat_outcome);
dat_outcome(isnan([dat_a';dat_b';dat_c';dat_d']))=NaN;
thresh_all(:,5)=[title1_low;num2cell(dat_outcome)];
thresh_all(:,10)=[[title1_low,' val'];num2cell([dat_a';dat_b';dat_c';dat_d'])];
thresh_all(:,15)=[[title1_low,'_fail_thresh_perc'];num2cell([graphm.thresh_fail.pos_perc{thresh_sel},graphm.thresh_fail.neg_perc{thresh_sel},...
    graphm.thresh_fail.pos_perc_test{thresh_sel},graphm.thresh_fail.neg_perc_test{thresh_sel}]')];


% % % participation coefficient
thresh_to_sel=20;thresh_pred=0.028056;true_pos = 'above';
thresh_sel=find(thresh_val == thresh_to_sel); % pruning level to select
dat_a=graphm.part_coef.part_pos_mean{thresh_sel};
dat_b=graphm.part_coef.part_neg_mean{thresh_sel};
dat_c=graphm.part_coef.part_pos_mean_test{thresh_sel};
dat_d=graphm.part_coef.part_neg_mean_test{thresh_sel};
title_1 = 'Participation coefficient';

% code below used for value export
if strcmp(true_pos,'above')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']>=thresh_pred;
else strcmp(true_pos,'below')
    dat_outcome=[dat_a';dat_b';dat_c';dat_d']<=thresh_pred;
end
title1_low =lower(title_1);
dat_outcome = double(dat_outcome);
dat_outcome(isnan([dat_a';dat_b';dat_c';dat_d']))=NaN;
thresh_all(:,6)=[title1_low;num2cell(dat_outcome)];
thresh_all(:,11)=[[title1_low,' val'];num2cell([dat_a';dat_b';dat_c';dat_d'])];

thresh_all(:,16)=[[title1_low,'_fail_thresh_perc'];num2cell([graphm.thresh_fail.pos_perc{thresh_sel},graphm.thresh_fail.neg_perc{thresh_sel},...
    graphm.thresh_fail.pos_perc_test{thresh_sel},graphm.thresh_fail.neg_perc_test{thresh_sel}]')];

% % write values to file
% xlswrite('D:\Google Drive\Arbeit\Lausanne\resting_state\Paper_Spectra\Stats\outcome_prediction_inc_spectra_inc_excluded_5.2-13.2.xlsx',thresh_all)

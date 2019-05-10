% Compute the frequency spectra and save them in respective folders
% copy (all) subjects to compute the freq spectra for from rs_psd_preprocessing 
%% 
% add required functions and toolboxes
addpath 'Path\to\power\spectra\scripts'
addpath 'Path\to\FieldTrip\'
ft_defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute frequency spectra %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_path = 'D:\Lausanne_analyses\RestingEEG\resting_data\'; % path in which patients are located
file_name = '\5s_segm_nothresh\data_interp.mat'; % path to move to within patient folder
out_file = 'data_freq_mtm.mat'; % output file name
elecs=p_layout('ladybird'); % type of cap being used
subj_miss_iter = 1;
for subj_iter = 1:numel(subjects)        
    display(['Loading data from: ',base_path,subjects{subj_iter},file_name]);
    load([base_path,subjects{subj_iter},file_name])
    % compute frequency spectra
    cfg                         = [];
    cfg.method                  = 'mtmfft';
    cfg.taper                   = 'dpss';
    cfg.tapsmofrq               = 1; % smoothing
    cfg.foi                     = 1:0.2:40;
%     cfg.keeptrials              = 'yes';
    frq                         = ft_freqanalysis(cfg,data);
    %
    mkdir([base_path,subjects{subj_iter},'5s_segm_nothresh\freq']);
    % Save results
    out_full = [base_path,subjects{subj_iter},'\5s_segm_nothresh\freq\',out_file];
    display(['Saving data to: ',out_full]);
%     save(out_file,'frq','lay');
    %
    keep subj_iter subjects base_path file_name keep elecs
end;

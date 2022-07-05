 %% Read BDF file
clear all
clc

% try
%      cd ('/Users/xxzhang/OneDrive - University College London/Ru_EEGSSVEP/TESTING/SCRIPTS_EEG_ANALYSIS')
% catch
%     cd ('C:\Users\lcharles\OneDrive - University College London\Experiments\Ru_EEGSSVEP\TESTING\SCRIPTS_EEG_ANALYSIS')%
%     %Lucies' laptop
% end


% add fieldtrip path
path_ft = '/Applications/MATLAB/toolbox/fieldtrip-20181231';
addpath(path_ft);
addpath('C:\Users\lcharles\Dropbox\WORK\MATLAB\fieldtrip-20191008')
addpath('D:\Dropbox\Dropbox\WORK\MATLAB\fieldtrip-20191008')


path_all = '/Users/xxzhang/OneDrive - University College London/Ru_EEGSSVEP/TESTING/';
addpath(genpath(path_all)); 
ft_defaults

% configuration

cfgPREPROC             = LZ_CONFIG_Analysis;


subjectCode            = cfgPREPROC.subjectCode;
Ntrials                = cfgPREPROC.Ntrials;
NArrows                = cfgPREPROC.nArrows;
Ndisplay               = cfgPREPROC.nDisplay;   

load('Biosemi64Layout.mat');

Nart = nan(length(subjectCode),1);
Ntri = nan(length(subjectCode),1);
trialrejection = nan(size((1:length(subjectCode)),2),1);



%% Reject trials with big variance ---------------------------------------------------

for subj = 1:length(subjectCode) % 1-12,15 no behav, 39 need redo
    disp(subj)
     % READ EEG FILE ---------------------------------------------------
    EEGfilesinfo = dir([cfgPREPROC.dir.data_EEGPreproc_ICAcleaned 'Subj_' num2str(subj,'%03.0f') '_Step3_ICACLEANED.mat']);
    load([cfgPREPROC.dir.data_EEGPreproc_ICAcleaned EEGfilesinfo.name]);
    
     % NEW BASELINE CORRECTION --------------------------------
     
         % -0.5  to 0 

         cfg.demean              = 'yes';
         cfg.baselinewindow      = [-0.500 0];
         dataNEW                 = ft_preprocessing(cfg,data);

         dataNEW.behav           = data.behav;
         
         % Correcting for sample info being from different
         % blocks (different lengths)
         for trial = 2:size(dataNEW.sampleinfo,1)
             
             if dataNEW.sampleinfo(trial,1) <= dataNEW.sampleinfo(trial-1,2)
                 dataNEW.sampleinfo(trial,1) = dataNEW.sampleinfo(trial,1) + dataNEW.sampleinfo(trial-1,2);
                 dataNEW.sampleinfo(trial,2) = dataNEW.sampleinfo(trial,2) + dataNEW.sampleinfo(trial-1,2);
             end
             
         end
                
         
    
%     cfg  = [];
%     cfg.baseline      = [-0.500 0];
%     dataNEW           = ft_timelockbaseline(cfg,data);
    

    % ARTIFACT DETECTION --------------------------------
    cfg = [];
    
        % function parameters

        cfg.continuous = 'no'; %z value can still be used, no matter continuous or not here
        %cfg.trl        = ?? %structure that defines the data segments of interest, see FT_DEFINETRIAL

        % channel selection, cutoff and padding parameters for artifact
        cfg.artfctdef.zvalue.channel      = 1:64; 
        cfg.artfctdef.zvalue.cutoff       = 20; % z-value threshold, stringent = between 4-20
        %cfg.artfctdef.zvalue.bpfilter     = 'yes';
        %cfg.artfctdef.zvalue.bpfreq       = [110 140]; %band pass frequency range (already previously bandpassed)
        cfg.artfctdef.zvalue.lpfilter      = 'yes';
        cfg.artfctdef.zvalue.lpfreq       = 48;
        %cfg.artfctdef.zvalue.bpfiltord    = 7; % related to smoothing
        %cfg.artfctdef.zvalue.bpfilttype   = 'but'; %digital filter type
        cfg.artfctdef.zvalue.hilbert      = 'yes';
        %cfg.artfctdef.zvalue.boxcar       = 0.2;
        cfg.artfctdef.zvalue.interactive = 'no'; % make the process interactive

        % what about:
            % cfg.artfctdef.zvalue.trlpadding = number in seconds
            % cfg.artfctdef.zvalue.fltpadding = number in seconds
            % cfg.artfctdef.zvalue.artpadding = number in seconds

         [cfg, artifact] = ft_artifact_zvalue(cfg,dataNEW); 
            % (not necessarily muscle artifact, so I changed artifact_muscle to artifact)
            % sensitive to EOG, muscle or jump artifacts; only works on continuously recorded data
            % scans data segments of interest for artifacts by means of thresholding the z-transformed value of the preprocessed raw data. 
   %%
    % ARTIFACT REJECTION --------------------------------
    cfg                                     = [];
    cfg.artfctdef.reject                    = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
    cfg.artfctdef.xxx.artifact              = artifact; 
    [data_no_artifacts, trialsremoved_art]  = ft_rejectartifact_NEW(cfg,dataNEW); 
        % output: cleaned data & automatically rejected trials
        
    trialskept_art = setdiff(1:length(dataNEW.trial),trialsremoved_art); % cleaned trials (remaining) after automatic rejection
     
    
    % FINAL INSPECTION --------------------------------
        % 
        % cfg                     = [];
        % cfg.layout              = layout;
        % cfg.method              = 'triangulation';
        % neighbours              = ft_prepare_neighbours(cfg);
        % 
        % cfg                     = [];
        % cfg.layout              = layout;
        % cfg.channel             = 1:64;
        % cfg.keepchannel         = 'repair';
        % cfg.neighbours          = neighbours;
        % cfg.keeptrial           = 'nan';
        % [data_trialcleaned, trialsremoved_vis] = ft_rejectvisual(cfg,data_no_artifacts); % why are there so many different names for the same data.
        % reject further trials by visual inspection
    
    % ALTERNATIVE PREPROCESSING STEP: no visual artifact detection
    trialsremoved_vis = [];
    
    
    % 5 - Saving removed trials
    % 1 for kept trials, 0 for trials to be removed
    trialsremove                        = [trialsremoved_art trialskept_art(trialsremoved_vis)]; % total trials to remove
    BADtrialsNEW                        = ones(length(data.trial),1); % create template structure with all trials
    BADtrialsNEW(trialsremove)          = 0; % all trials with total bad trials removed
    BADtrials_artNEW                    = ones(length(data.trial),1);  
    BADtrials_artNEW(trialsremoved_art) = 0; % all trials with muscle artifact bad trials removed
    BADtrials_visNEW                    = ones(length(data.trial),1);
    BADtrials_visNEW(trialskept_art(trialsremoved_vis)) = 0; % all trials with visual bad trials removed
    
    Nart(subj) = size(artifact,1);
    Ntri(subj) = size(trialsremoved_art,2);
    
    %Save in the original file
    filename = [cfgPREPROC.dir.data_EEGPreproc_FINAL 'Subj_' num2str(subj,'%03.0f') '_FINAL_Z20.mat'];
    save(filename,'dataNEW','-v7.3') % why -v7.3?
    save(filename,'-append','BADtrialsNEW','BADtrials_artNEW','BADtrials_visNEW');

    clear data_trialcleaned data_no_artifacts data    
    
    % calculate and save trial rejection rate
    trialrejection(subj) = Ntri(subj)./trial;

end

% ave_trialrejection = nanmean(trialrejection);
% filename = [cfgPREPROC.dir.data_EEGPreproc_FINAL 'ALLSUBJ_autorejrate.mat'];
% save(filename,'trialrejection','ave_trialrejection');

%% Checking artifacts detected & trials removed
% 
% for subj = 1:5 %1:length(subjectCode)
%     disp(subj)
%      % READ EEG FILE ---------------------------------------------------
%     EEGfile = dir([cfgPREPROC.dir.data_EEGPreproc_FINAL,'Subj_', num2str(subj,'%03.0f'), '_FINAL_Z20.mat']);
%     load([cfgPREPROC.dir.data_EEGPreproc_FINAL EEGfile.name]);
%     
%     Nart(subj) = size(artifact,1);
%     Ntri(subj) = size(trialsremoved_art,2);
% 
% end
%% DRAFT
% 
%       
%       cfg=[];
%       cfg.artfctdef.muscle.channel     = 1:64;
%       cfg.continuous='no';
%       cfg.artfctdef.muscle.bpfilter    = 'yes';
%       cfg.artfctdef.muscle.bpfreq      = [150 180];
%       cfg.artfctdef.muscle.bpfiltord   = 7;
%       cfg.artfctdef.muscle.bpfilttype  = 'but';
%       cfg.artfctdef.muscle.hilbert     = 'yes';
% %       cfg.artfctdef.muscle.trlpadding  = 0.1;
% %       cfg.artfctdef.muscle.fltpadding  = 0.1;
% %       cfg.artfctdef.muscle.artpadding  = 0.1;
% 
%         [cfg, artifactm] = ft_artifact_muscle(cfg, data);



%  % IF no extra trial deleted - Saving removed trials
%      trialsremove                        = [trialsremoved_art];
%      BADtrialsNEW                           = ones(length(data.trial),1);
%      BADtrialsNEW(trialsremove)             = 0;
%      BADtrials_artNEW                       = ones(length(data.trial),1);
%      BADtrials_artNEW(trialsremoved_art)  	= 0;
%      
%      
%      sum(BADtrials_artNEW==0)
%      sum(BADtrials_visNEW==0)

%% This is Lemona's script for:
    % Frequency analysis
    % trial-by-trial SSVEP & behavioural data saved

%% SET UP

clear
clc
try
    cd ('/Users/xxzhang/OneDrive - University College London/Ru_EEGSSVEP/TESTING/SCRIPTS_EEG_ANALYSIS')
catch
%     cd ('C:\Users\lcharles\OneDrive - University College London\Experiments\Ru_EEGSSVEP\TESTING\SCRIPTS_EEG_ANALYSIS')
end


%add fieldtrip path
path_ft = '/Applications/MATLAB/toolbox/fieldtrip-20181231';
addpath(path_ft);
% addpath('D:\Dropbox\Dropbox\WORK\MATLAB\fieldtrip-2011')
addpath('D:\Dropbox\Dropbox\WORK\MATLAB\fieldtrip-20191008')
addpath('C:\Users\lcharles\Dropbox\WORK\MATLAB\fieldtrip-20191008')

path_all = '/Users/xxzhang/OneDrive - University College London/Ru_EEGSSVEP/TESTING/';
addpath(genpath(path_all)); % cannot find LZ_CONFIG_Analysis if set path to Step 4...
ft_defaults

addpath('../')

cfgPREPROC             = LZ_CONFIG_Analysis;

Ntrials                = cfgPREPROC.Ntrials;
NArrows                = cfgPREPROC.nArrows;
Ndisplay               = cfgPREPROC.nDisplay;


load('Biosemi64Layout.mat');

%% DATA SETUP

subjlist            = cfgPREPROC.subjectCode; 
%subjlist = setdiff(subjlist,[1 2 7 8 10 18 35 37 38 40 50 53]); 
    % remove bad subjects identified by eyeballing averaged ERPs
occ_channels    = 56:64;
front_channels  = 1:17;
frequencybins   = [5 10];
    
GOODSUBJECTS = subjlist;

%% Load data & time frequency analysis (save trial-by-trial frequency & behavioural info)    
% data structure: number of subjects x number of channels x 2 frequencies
MEAN_freq_OCC = nan(length(subjlist),size(occ_channels,2),2);
STD_freq_OCC = nan(length(subjlist),size(occ_channels,2),2);
MEAN_freq_FRONT = nan(length(subjlist),size(front_channels,2),2);
STD_freq_FRONT = nan(length(subjlist),size(front_channels,2),2);


for subj = GOODSUBJECTS
    disp(subj)

    EEGfile = dir([cfgPREPROC.dir.data_EEGPreproc_FINAL,'Subj_', num2str(subj,'%03.0f'), '_FINAL_Z20.mat']);
%     filename = strcat('Subj_', num2str(subj,'%03.0f'),'_FINAL_Z20.mat');
    %EEG = pop_loadset('filename',filename,'filepath', path_data);
    load([cfgPREPROC.dir.data_EEGPreproc_FINAL EEGfile.name]);
    % load("Subj_001_FINAL.mat");
%     load(filename);

    behav            = dataNEW.behav;
    behav = rmfield(behav,{'angles'}); % remove useless field from the behavioural data structure
    behav = rmfield(behav,{'respgoldenArrow'});
    filterkepttrials = BADtrialsNEW;

    %4/ analysis 4 :Stimulus locked TFR : specific freq
    % - 5 Hz / 10 Hz /
    % - Occ channels / frontal channels/

    for thisfreq = 1:2
        filtertrials = behav.GoldenTrial == 0 & behav.resp > 0;
        
        if subj == 21 
            filtertrials(end - 1:end) = []; % remove last 2 trials
        elseif subj == 46
            filtertrials(1:2) = []; % remove first 2 trials
        end
        
        %trials analysed: 
        % not golden arrow attention check trials 
        % behavioural response is not 0
        %behavioural response coding:0-no response, 1-left, 2-right

        %OCCIPITAL CHANNELS

        cfg = [];
        cfg.method      = 'mtmfft';
        %mtmfft performs frequency analysis on any time series trial data using a conventional
        % single taper (e.g. Hanning) or using the multiple tapers based on discrete prolate
        % spheroidal sequences (DPSS), also known as the Slepian sequence.
        cfg.pad         = 'nextpow2'; %length in seconds to which the data can be padded out. The
        %padding will determine your spectral resolution.
        cfg.foi         = frequencybins(thisfreq); %vector 1 x numfoi, frequencies of interest
        cfg.channel     = occ_channels;
        cfg.taper       = 'hanning';
        cfg.output      = 'powandcsd';
        cfg.keeptrials  = 'yes';
        cfg.toi         = 0.25:5.45;  %vector 1 x numtoi, the times on which the analysis windows should be centered (in seconds)
        cfg.trials      = find(filtertrials' & filterkepttrials);
        freq_analysis   = ft_freqanalysis(cfg, dataNEW);

        
        
        % Apply baseline correction and save in structure
        thismean    = mean(freq_analysis.powspctrm,1);
        thissdt     = std(freq_analysis.powspctrm);
        dataz       = (freq_analysis.powspctrm - thismean )./ thissdt;       
        SSVEP(subj).freq_trialbytrial{thisfreq,1} = dataz; % remaining trials x n channels (per location) data points per pp per frequency
        

      
        %FRONTAL CHANNELS
        cfg = [];
        cfg.method      = 'mtmfft';
        cfg.pad         = 'nextpow2';
        cfg.foi         = frequencybins(thisfreq);
        cfg.channel     = front_channels;
        cfg.taper       = 'hanning';
        cfg.output      = 'powandcsd';
        cfg.toi         = 0.25:5.45;
        cfg.keeptrials  = 'yes';
        cfg.trials      = find(filtertrials' & filterkepttrials);
        freq_analysis   = ft_freqanalysis(cfg, dataNEW);
        
        
        % Apply baseline correction and save in structure
        thismean    = mean(freq_analysis.powspctrm,1);
        thissdt     = std(freq_analysis.powspctrm);
        dataz       = (freq_analysis.powspctrm - thismean )./ thissdt;       
        SSVEP(subj).freq_trialbytrial{thisfreq,2} = dataz; % 336 trials x n channels (per location) data points per pp per frequency
         
        
        % save trial by trial frequency & behavioural info in table format
        tmptable = struct2table(behav,'AsArray',1); % create a temporary table that converts the structure into table (more easily manipulable)
        A = []; % this should have 384/336 trials?  
        for thiscol = 1:size(tmptable,2) % for each behavioural property (column)
            thisppty = tmptable{1,thiscol}; % for this property (column), define value of all trials
            A = [A thisppty(cfg.trials )']; % put value of all trials into SSVEP
        end
        SSVEP(subj).behav = array2table(A,'VariableNames',fieldnames(behav));
        
        clear freq_analysis
    end
    

end

%% Check overall trial rejection rate

% Load trial-by-trial data from FINAL Script 2
% filename = 'ALLSUBJ_freqbehav_trialbytrial.mat';
% load([cfgPREPROC.dir.data_EEG_analysis filename])

% Summarise trials removed and calculate average trial rejection rate
trialrejection = nan(62,1);
% sum(BADtrialsNEW) from this script.

for subj = subjlist
trialrejection(subj) = 1 - (size(SSVEP(subj).freq_trialbytrial{1,1},1) ./ ... 
size(filterkepttrials,1));
end

ave_trialrejection = nanmean(trialrejection);

%% Save SSVEP analysed data and trial rejection info

save([cfgPREPROC.dir.data_EEG_analysis 'ALLSUBJ_freqbehav_trialbytrial_v2.mat'],'SSVEP','trialrejection','ave_trialrejection'); 


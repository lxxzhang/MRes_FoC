%% This is Lemona's frequency analysis script.

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

%% LOAD DATA
subj_EEG     = setdiff(cfgPREPROC.subjectCode,[3 5 11 13 14 16 17 19 20 23 25 27 28 29 32 33 34 36 39 45 47 48 51 56 58 59 61 62]); 
subjlist        = cfgPREPROC.subjectCode;
occ_channels    = 56:64;
front_channels  = 1:17;
frequencybins   = [5 10];

% LOAD DATA (SSVEP: trial-by-trial frequency & behavioural data)
load([cfgPREPROC.dir.data_EEG_analysis 'ALLSUBJ_freqbehav_trialbytrial_v2.mat']); 

%% ANALYSES BY CONDITION X CONGRUENCE

 for subj = subjlist 
     for thisregion = 1:2
        for thisfreq = 1:2
            for instr = 1:3
                for instronset = 1:2
                    %% TRIAL SELECTION
                    % select relevant trials by condition (corresponding to current condition
                        % combination & non-golden)
                        filter_cond = SSVEP(subj).behav.responseMode == instr & ...
                        SSVEP(subj).behav.instrOnset == instronset & ...
                        SSVEP(subj).behav.GoldenTrial == 0 & SSVEP(subj).behav.resp > 0;
  
                    % select congruent/incongruent trials
                    for respcong = 1:2
                        % select congruent trials
                        if respcong == 1 
                            filter_cong = SSVEP(subj).behav.orientation == SSVEP(subj).behav.resp;
                            
                        % select incongruent trials
                        else
                            filter_cong = SSVEP(subj).behav.orientation ~= SSVEP(subj).behav.resp;
                        end
                        % congruence filter on top of basic condition filter
                        filter_cong = filter_cond & filter_cong;
                        
                        % check trials assigned correctly
%                         if thisfreq == 1 & thisregion == 1
%                             disp(['Subj' num2str(subj) ' Instr ' num2str(instr) ' Onset ' num2str(instronset) ' acc is ' num2str(sum(filter_cong)/sum(filter_cond))])
%                         end
                        
                        % trial limit filter (only by condition)
                        Ntriallimit = 20;
                        if sum(filter_cond) < Ntriallimit
                            data_cong(subj,instr,instronset,respcong,thisfreq,thisregion) = 0; % discard subjects with too few relevant trials
                        
                       else
                            % analysis by condition x congruence (average across trials)
                            data_cong(subj,instr,instronset,respcong,thisfreq,thisregion)  = mean(nanmean(SSVEP(subj).freq_trialbytrial{thisfreq,thisregion}(filter_cong,:),2));
                            
                            %code accuracy data based on instruction type
                            if instr == 1 || instr == 3
                                data_acc(subj,instr,instronset,respcong,thisfreq,thisregion) = data_cong(subj,instr,instronset,respcong,thisfreq,thisregion);
                            else
                                if respcong == 1
                                    data_acc(subj,instr,instronset,2,thisfreq,thisregion) = data_cong(subj,instr,instronset,respcong,thisfreq,thisregion);
                                else
                                    data_acc(subj,instr,instronset,1,thisfreq,thisregion) = data_cong(subj,instr,instronset,respcong,thisfreq,thisregion);
                                end
                            end
                        end
                    end
                end
            end
        end
     end
 end
            
save([cfgPREPROC.dir.data_EEG_analysis 'SSVEP_cong_N34.mat'],'data_cong'); 
            
%% Plot mean SSVEP (violin) by condition x congruence
% subj,instr,instronset,thisfreq,thisregion

%% 5Hz
figure;
thisfreq = 1;
plotfreq = 5;

% data_cong: subj,instr,instronset,respcong,freq,region

for thisregion = 1:2
     if thisregion == 1
            nameregion = ' Occipital ';
        elseif thisregion == 2
            nameregion = ' Frontal ';
     end    
     % Remove any participant with any zeros value in any condition
    GOODSUBJECTS = subjlist;
    for instr = 1:3
        for instronset = 1:2
            for respcong = 1:2
                GOODSUBJforthiscond = any(data_cong(GOODSUBJECTS,instr,instronset,respcong,thisfreq,thisregion),2);
                GOODSUBJECTS = GOODSUBJECTS(GOODSUBJforthiscond);
            end
        end
    end
    for instr = 1:3
        subplot(2,3,instr + (thisregion-1)*3)
        
        % each pp's data for each condition
        % condition x congruent x before/after
        y = squeeze(data_cong(GOODSUBJECTS,instr,:,1,thisfreq,thisregion));
        
        % condition x incongruent x before/after
        z = squeeze(data_cong(GOODSUBJECTS,instr,:,2,thisfreq,thisregion));
        if instr == 1
            colour_cong = [0.73 1 0.31; 0.54 0.9 0]; % cong before, cong after
            colour_incong = [0.51 0.7 0.44; 0.31 0.6 0.2]; % incong before, incong after
            nameinstr = 'Adhere ';
        elseif instr == 2
            colour_cong = [1 0.55 0.26; 0.85 0.36 0]; % cong before, cong after
            colour_incong = [0.93 0.3 0.3; 0.78 0.15 0.15]; % incong before, incong after
            nameinstr = 'Oppose ';
        else
            colour_cong = [0.3 0.91 1; 0 0.66 0.76]; % cong before, cong after
            colour_incong = [0.29 0.67 1; 0 0.39 0.73]; % incong before, incong after
            nameinstr = 'Detach ';
        end
        
        al_goodplot(y,1:2,0.5,colour_cong,'left',[],[]); % Congruent
        al_goodplot(z,1:2,0.5,colour_incong,'right',[],[]); % Incongruent
        grid on
        ylim([-0.5 0.7])
        xlim([0.3 2.7])
        xticks([1 2])
        xticklabels({'Before', 'After'})
        
        thistitle = sprintf(strcat(nameinstr, nameregion));
        title(thistitle);
        
        % T-TESTS: for before or after timing, t-test between congruent and incongruent SSVEPs
                for instronset = 1:2
                    thisonset = [' Before';'After  '];
                    [res_str,p,d,resshort_str] = disp_tstat5_CohenD_APA(y(:,instronset), ...
                        z(:,instronset),'Congruent','Incongruent', 1,'notail',thisonset(instronset,:),thistitle);
                    if p <= 0.05
                        sigstar([instronset-0.25 instronset+0.25],0.05)
        
                        % 1 - ispaired
                        %Freq Power = variable name (?)
                    end
                end
    end
end


%% 10Hz
figure;
thisfreq = 2;
plotfreq = 10;

% data_cong: subj,instr,instronset,respcong,freq,region

for thisregion = 1:2
     if thisregion == 1
            nameregion = ' Occipital ';
        elseif thisregion == 2
            nameregion = ' Frontal ';
     end    
     % Remove any participant with any zeros value in any condition
    GOODSUBJECTS = subjlist;
    for instr = 1:3
        for instronset = 1:2
            for respcong = 1:2
                GOODSUBJforthiscond = any(data_cong(GOODSUBJECTS,instr,instronset,respcong,thisfreq,thisregion),2);
                GOODSUBJECTS = GOODSUBJECTS(GOODSUBJforthiscond);
            end
        end
    end
    for instr = 1:3
        subplot(2,3,instr + (thisregion-1)*3)
        
        % each pp's data for each condition
        % condition x congruent x before/after
        y = squeeze(data_cong(GOODSUBJECTS,instr,:,1,thisfreq,thisregion));
        
        % condition x incongruent x before/after
        z = squeeze(data_cong(GOODSUBJECTS,instr,:,2,thisfreq,thisregion));
        if instr == 1
            colour_cong = [0.73 1 0.31; 0.54 0.9 0]; % cong before, cong after
            colour_incong = [0.51 0.7 0.44; 0.31 0.6 0.2]; % incong before, incong after
            nameinstr = 'Adhere ';
        elseif instr == 2
            colour_cong = [1 0.55 0.26; 0.85 0.36 0]; % cong before, cong after
            colour_incong = [0.93 0.3 0.3; 0.78 0.15 0.15]; % incong before, incong after
            nameinstr = 'Oppose ';
        else
            colour_cong = [0.3 0.91 1; 0 0.66 0.76]; % cong before, cong after
            colour_incong = [0.29 0.67 1; 0 0.39 0.73]; % incong before, incong after
            nameinstr = 'Detach ';
        end
        
        al_goodplot(y,1:2,0.5,colour_cong,'left',[],[]); % Congruent
        al_goodplot(z,1:2,0.5,colour_incong,'right',[],[]); % Incongruent
        grid on
        ylim([-0.5 0.7])
        xlim([0.3 2.7])
        xticks([1 2])
        xticklabels({'Before', 'After'})
        
        thistitle = sprintf(strcat(nameinstr, nameregion));
        title(thistitle);
        
        % T-TESTS: for before or after timing, t-test between congruent and incongruent SSVEPs
                for instronset = 1:2
                    thisonset = [' Before';'After  '];
                    [res_str,p,d,resshort_str] = disp_tstat5_CohenD_APA(y(:,instronset), ...
                        z(:,instronset),'Congruent','Incongruent', 1,'notail',thisonset(instronset,:),thistitle);
                    if p <= 0.05
                        sigstar([instronset-0.25 instronset+0.25],0.05)
        
                        % 1 - ispaired
                        %Freq Power = variable name (?)
                    end
                end
    end
end

%% ANOVAs
    % data_cong:subj,instr,instronset,respcong,freq,region
    % ADHERE and OPPOSE: instruction, instruction onset, response accuracy

    
disp(['ANOVA (congruence) - SSVEP for ADHERE and OPPOSE conditions for at 5Hz at OCCIPITAL channels']);
d = data_cong(GOODSUBJECTS,1:2,:,:,1,1);
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Congruence'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA (accuracy) - SSVEP for ADHERE and OPPOSE conditions at 5Hz at OCCIPITAL channels']);
d = data_acc(GOODSUBJECTS,1:2,:,:,1,1);
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Accuracy'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA (congruence) - SSVEP for ADHERE and OPPOSE conditions at 5Hz at OCCIPITAL channels BEFORE']);
d = squeeze(data_cong(GOODSUBJECTS,1:2,1,:,1,1));
instr_anova = simple_mixed_anova(d, [],{'Instr','Congruence'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA (congruence) - SSVEP for ADHERE and OPPOSE conditions at 5Hz at OCCIPITAL channels AFTER']);
d = squeeze(data_cong(GOODSUBJECTS,1:2,2,:,1,1));
instr_anova = simple_mixed_anova(d, [],{'Instr','Congruence'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA (congruence) - SSVEP for ADHERE and OPPOSE conditions at 10Hz at OCCIPITAL channels']);
d = data_cong(GOODSUBJECTS,1:2,:,:,2,1);
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Congruence'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA (accuracy) - SSVEP for ADHERE and OPPOSE conditions for all channels at 10Hz at OCCIPITAL channels']);
d = data_acc(GOODSUBJECTS,1:2,:,:,2,1);
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Accuracy'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

% 10Hz fr

disp(['ANOVA (congruence) - SSVEP for ADHERE and OPPOSE conditions at 5Hz at FRONTAL channels']);
d = data_cong(GOODSUBJECTS,1:2,:,:,1,2);
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Congruence'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

% disp(['ANOVA (accuracy) - SSVEP for ADHERE and OPPOSE conditions for all channels at 5Hz at FRONTAL channels']);
% d = data_acc(GOODSUBJECTS,1:2,:,:,1,2);
% instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Accuracy'}); % labels for each variable
% instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA (congruence) - SSVEP for ADHERE and OPPOSE conditions at 5Hz at ALL channels']);
d = squeeze(data_cong(GOODSUBJECTS,1:2,:,:,1,1:2));
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Congruence','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA (accuracy) - SSVEP for ADHERE and OPPOSE conditions at 5Hz at ALL channels']);
d = squeeze(data_acc(GOODSUBJECTS,1:2,:,:,1,1:2));
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Accuracy','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

% DETACH only 
disp(['ANOVA - SSVEP for DETACH condition at 5Hz at OCCIPITAL channels']);
d = squeeze(data_cong(GOODSUBJECTS,3,:,:,1,1));
instr_anova = simple_mixed_anova(d, [],{'InstrOnset','Congruence'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - SSVEP for DETACH condition at 5Hz at ALL channels']);
d = squeeze(data_cong(GOODSUBJECTS,3,:,:,1,:));
instr_anova = simple_mixed_anova(d, [],{'InstrOnse. bv    t','Congruence','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - SSVEP for DETACH condition at 10Hz at ALL channels']);
d = squeeze(data_cong(GOODSUBJECTS,3,:,:,2,:));
instr_anova = simple_mixed_anova(d, [],{'InstrOnset','Congruence','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

%% pasted from tbt analysis script
% disp(['ANOVA SSVEP FOR ADHERE AND OPPOSE conds FOR ALL channels AT 5HZ']);
% d = squeeze(data_cond(GOODSUBJECTS,1:2,:,1,:));
% instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Region'}); % labels for each variable
% instr_anova = partial_eta_squared(instr_anova)
% 
% 
% disp(['ANOVA SSVEP FOR DETACH conds FOR ALL channels']);
% d = squeeze(data_cond(GOODSUBJECTS,3,:,:,:));
% instr_anova = simple_mixed_anova(d, [],{'InstrOnset','Freq','Region'}); % labels for each variable
% instr_anova = partial_eta_squared(instr_anova)
% 
% 
% disp(['ANOVA SSVEP FOR DETACH conds FOR ALL channels AT 5HZ']);
% d = squeeze(squeeze(data_cond(GOODSUBJECTS,3,:,1,:)));
% instr_anova = simple_mixed_anova(d, [],{'InstrOnset','Region'}); % labels for each variable
% instr_anova = partial_eta_squared(instr_anova)
% 
% 
% disp(['ANOVA SSVEP FOR ADHERE AND DETACH conds FOR ALL channels']);
% d = data_cond(GOODSUBJECTS,[1 3],:,:,:);
% instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Freq','Region'}); % labels for each variable
% instr_anova = partial_eta_squared(instr_anova)
% 
% 
% disp(['ANOVA SSVEP FOR ADHERE AND DETACH conds FOR ALL channels AT 5HZ']);
% d = squeeze(data_cond(GOODSUBJECTS,[1 3],:,1,:));
% instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Region'}); % labels for each variable
% instr_anova = partial_eta_squared(instr_anova)

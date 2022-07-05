%% This is Lemona's frequency analysis script by congruence.

%% SET UP

% clear
% clc
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

subjlist            = cfgPREPROC.subjectCode;

occ_channels    = 56:64;
front_channels  = 1:17;
frequencybins   = [5 10];

% LOAD DATA (SSVEP: trial-by-trial frequency & behavioural data)
load([cfgPREPROC.dir.data_EEG_analysis 'ALLSUBJ_freqbehav_trialbytrial_v2.mat']); 


%% ANALYSES BY CONDITION

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
                    
                    % analysis by condition (redo previous sanity check) 
                    Ntriallimit = 20;
                    if sum(filter_cond) < Ntriallimit
                        data_cond(subj,instr,instronset,thisfreq,thisregion) = 0; % discard subjects with too few relevant trials
                    else
                        data_cond(subj,instr,instronset,thisfreq,thisregion)  = mean(mean(SSVEP(subj).freq_trialbytrial{thisfreq,thisregion}(filter_cond,:),2)); % average across channels across trials
                    end
                end
            end
        end
     end
 end
            
            
%% Plot mean SSVEP (violin) by condition
% subj,instr,instronset,thisfreq,thisregion

figure;
frequencieslist = [5 10];
% data_cond: subj,instr,instronset,freq,region

for thisregion = 1:2
    if thisregion == 1
        nameregion = ' Occipital ';
    elseif thisregion == 2
        nameregion = ' Frontal ';
    end  
    for thisfreq = 1:2
  
        plotfreq = frequencieslist(thisfreq);
    
        subplot(2,2,thisfreq + (thisregion-1)*2)
        % FRONT 5Hz, FRONT 10Hz, OCC 5Hz, OCC 10Hz

        % Remove any participant with any zeros value in any condition
        GOODSUBJECTS = subjlist;
            for instr = 1:3
                for instronset = 1:2
                     GOODSUBJforthiscond = any(data_cond(GOODSUBJECTS,instr,instronset,thisfreq,thisregion),2);
                     GOODSUBJECTS = GOODSUBJECTS(GOODSUBJforthiscond);
                end
            end

        % each pp's data for each condition
        y = data_cond(GOODSUBJECTS,:,1,thisfreq,thisregion); % before
        z = data_cond(GOODSUBJECTS,:,2,thisfreq,thisregion); % after
        
        % (colour name switched)
        al_goodplot(y,1:3,0.5,cfgPREPROC.colorinstrAFTER,'left',[5 5 5],[]); % Before
        al_goodplot(z,1:3,0.5,cfgPREPROC.colorinstrBEFORE,'right',[5 5 5],[]); % After

        ylim([-0.5 0.5])
        xlim([0.5 3.5])
        xticks([1 2 3])
        xticklabels({'Adhere', 'Oppose', 'Detach'})
        colorinstrBEFORE = cfgPREPROC.colorinstrBEFORE;
        colorinstrAFTER = cfgPREPROC.colorinstrAFTER;
       
        thistitle = sprintf(strcat(num2str(plotfreq),' Hz',nameregion));
        title(thistitle);

        %% T-TESTS: for each instruction type, t-test between overall instruction timing
        for instr = 1:3
            thiscond = [' Adhere';'Oppose ';'Detach '];
            [res_str,p,d,resshort_str] = disp_tstat5_CohenD_APA(y(:,instr), ...
                z(:,instr),'BEFORE','AFTER', 1,'notail',thiscond(instr,:),thistitle);
            if p <= 0.05
                sigstar([instr-0.25 instr+0.25],0.05)

                % 1 - ispaired
                %Freq Power = variable name (?)
            end
        end
    end
end

%% ANOVAs
    % data_cond: subj,instr,instronset,freq,region
    
disp(['ANOVA SSVEP FOR ADHERE AND OPPOSE conds FOR ALL channels']);
d = data_cond(GOODSUBJECTS,1:2,:,:,:);
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Freq','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA SSVEP FOR ADHERE AND OPPOSE conds FOR ALL channels AT 5HZ']);
d = squeeze(data_cond(GOODSUBJECTS,1:2,:,1,:));
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)


disp(['ANOVA SSVEP FOR DETACH conds FOR ALL channels']);
d = squeeze(data_cond(GOODSUBJECTS,3,:,:,:));
instr_anova = simple_mixed_anova(d, [],{'InstrOnset','Freq','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)


disp(['ANOVA SSVEP FOR DETACH conds FOR ALL channels AT 5HZ']);
d = squeeze(squeeze(data_cond(GOODSUBJECTS,3,:,1,:)));
instr_anova = simple_mixed_anova(d, [],{'InstrOnset','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)


disp(['ANOVA SSVEP FOR ADHERE AND DETACH conds FOR ALL channels']);
d = data_cond(GOODSUBJECTS,[1 3],:,:,:);
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Freq','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)


disp(['ANOVA SSVEP FOR ADHERE AND DETACH conds FOR ALL channels AT 5HZ']);
d = squeeze(data_cond(GOODSUBJECTS,[1 3],:,1,:));
instr_anova = simple_mixed_anova(d, [],{'Instr','InstrOnset','Region'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

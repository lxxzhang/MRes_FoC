%% This is Lemona's script for reverse correlations and plotting.

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

cfgPREPROC = LZ_CONFIG_Analysis;

subj_EEG   = setdiff(cfgPREPROC.subjectCode,[3 5 11 13 14 16 17 19 20 23 25 27 28 29 32 33 34 36 39 45 47 48 51 56 58 59 61 62]); 
subjlist   = subj_EEG;
Ntrials    = cfgPREPROC.Ntrials;
NArrows    = cfgPREPROC.nArrows;
Ndisplay   = cfgPREPROC.nDisplay;


load('Biosemi64Layout.mat');

%% LOAD DATA

occ_channels    = 56:64;
front_channels  = 1:17;
frequencybins   = [5 10];

% LOAD DATA (SSVEP: trial-by-trial frequency & behavioural data)
load([cfgPREPROC.dir.data_EEG_analysis 'ALLSUBJ_freqbehav_trialbytrial_v2.mat']); 
    
%% response congruence rate X SSVEP amplitude (median split)
for subj = subjlist
    for thisregion = 1 % occipital
        for thisfreq = 2 % 10Hz
            for instronset = 1:2
                for instr = 1:3
                    
                   filter_cond= SSVEP(subj).behav.responseMode == instr & ...
                        SSVEP(subj).behav.instrOnset == instronset & ...
                        SSVEP(subj).behav.GoldenTrial == 0 & SSVEP(subj).behav.resp > 0;
                    
                    % compute the median SSVEP for that condition
                    %MedianSSVEP = median(mean(SSVEP(subj).freq_trialbytrial{thisfreq,thisregion}(filter_cond,:),2));
                    MedianSSVEP = median(mean(SSVEP(subj).freq_trialbytrial{thisfreq,thisregion}(SSVEP(subj).behav.GoldenTrial == 0 & SSVEP(subj).behav.resp > 0,:),2));
                        % Grand median for all conditions
                    
                    for SSVEPamp =1:2
                        if SSVEPamp == 1
                            filterSSVEPamp = mean(SSVEP(subj).freq_trialbytrial{thisfreq,thisregion}(filter_cond,:),2) <= MedianSSVEP;                        
                            % average across channels

%                             SSVEPampLO = filterSSVEPamp .* SSVEP(subj).freq_trialbytrial{thisfreq,thisregion};
%                             SSVEPampLO(SSVEPampLO==0)=NaN;
%                             aveSSVEPampLO = mean(nanmean(SSVEPampLO));
                        elseif SSVEPamp == 2
                            filterSSVEPamp = mean(SSVEP(subj).freq_trialbytrial{thisfreq,thisregion}(filter_cond,:),2) > MedianSSVEP;
%                             SSVEPampHI = filterSSVEPamp .* SSVEP(subj).freq_trialbytrial{thisfreq,thisregion};
%                             SSVEPampHI(SSVEPampHI==0)=NaN;
%                             aveSSVEPampHI = mean(nanmean(SSVEPampHI));
                            
                        end
                  
                        cong_trials = size(find(SSVEP(subj).behav.orientation(filterSSVEPamp,:) == SSVEP(subj).behav.resp(filterSSVEPamp,:)),1);
                            % for the participants with higher/lower SSVEP
                            % amplitude, find the number of congruent
                            % trials within this condition combination.
                        cong_rate(subj,SSVEPamp,instr,instronset) = cong_trials ./ size(find(filterSSVEPamp),1);
                            % last remaining subject number is 60
                        cong_rate(cong_rate == 0) = nan;
                        
                end
            end
        end
    end
end
end

%% Plot

f1 = figure('Name','Response congruence rate X SSVEP amplitude group','units','normalized','outerposition',[0 0 0.55 0.90]);

cong_rate_ave = squeeze(nanmean(cong_rate,1));
    %SSVEPamp,instr,instronset
for instr = 1:3
        disp(instr)
        set(gcf,'color','white');
        subplot(1,3,instr);
        bar([0.5 1],cong_rate_ave(1:2,instr,1),'FaceColor',cfgPREPROC.colorinstrAFTER(instr,:),'EdgeColor',[1 1 1],'BarWidth',0.8);hold on;
        bar([2 2.5],cong_rate_ave(1:2,instr,2),'FaceColor',cfgPREPROC.colorinstrBEFORE(instr,:),'EdgeColor',[1 1 1],'BarWidth',0.8);hold on;

        %errorbar([1:2],Accuracy_condave(instr,:),Accuracy_condSE(instr,:),'LineWidth',3,'color',[0.8 0.8 0.8],'LineStyle','none');
     
        
        xlabel('BEFORE       AFTER','FontWeight','bold','FontSize',15);

        xlim([0 3]); 
        xticks([0.5 2]);
        set(gca,'xtick',[0.75 2.25],'XtickLabel','Low  High','FontSize',12); 
        

        ylabel('% Congruent Responses','FontWeight','bold','FontSize',16);
        
        ylim([0 1]);
        set(refline(0,0.5),'LineStyle','--','color','k');
        title(cfgPREPROC.InstrLabels{instr},'FontSize',18,'Color',cfgPREPROC.colours(instr,:))
        box off;
end

%% ANOVAs
% all at 10Hz and occipital channel
% cong_rate(subj,SSVEPamp,instr,instronset)

% ADHERE x OPPOSE
disp(['ANOVA - congruence x SSVEP amplitude for ADHERE and OPPOSE conditions']);
d = cong_rate(subjlist,:,1:2,:); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp','Instruction','InstrOnset'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for ADHERE and OPPOSE conditions for BEFORE instruction onset']);
d = squeeze(cong_rate(subjlist,:,1:2,1)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp','Instruction'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for ADHERE and OPPOSE conditions for AFTER instruction onset']);
d = squeeze(cong_rate(subjlist,:,1:2,2)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp','Instruction'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

% DETACH
disp(['ANOVA - congruence x SSVEP amplitude for DETACH condition for all instruction onset']);
d = squeeze(cong_rate(subjlist,:,3,:)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp','InstrOnset'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for DETACH condition for BEFORE instruction onset']);
d = squeeze(cong_rate(subjlist,:,3,1)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp'}); 
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for DETACH condition for AFTER instruction onset']);
d = squeeze(cong_rate(subjlist,:,3,2)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

% ADHERE only
disp(['ANOVA - congruence x SSVEP amplitude for ADHERE condition for all instruction onset']);
d = squeeze(cong_rate(subjlist,:,1,:)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp','InstrOnset'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for ADHERE condition for BEFORE instruction onset']);
d = squeeze(cong_rate(subjlist,:,1,1)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp'}); 
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for ADHERE condition for AFTER instruction onset']);
d = squeeze(cong_rate(subjlist,:,1,2)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

% OPPOSE only
disp(['ANOVA - congruence x SSVEP amplitude for OPPOSE condition for all instruction onset']);
d = squeeze(cong_rate(subjlist,:,2,:)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp','InstrOnset'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for OPPOSE condition for BEFORE instruction onset']);
d = squeeze(cong_rate(subjlist,:,2,1)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp'}); 
instr_anova = partial_eta_squared(instr_anova)

disp(['ANOVA - congruence x SSVEP amplitude for OPPOSE condition for AFTER instruction onset']);
d = squeeze(cong_rate(subjlist,:,2,2)); 
instr_anova = simple_mixed_anova(d, [],{'SSVEPamp'}); % labels for each variable
instr_anova = partial_eta_squared(instr_anova)


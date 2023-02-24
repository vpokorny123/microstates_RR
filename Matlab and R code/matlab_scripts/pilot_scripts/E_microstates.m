addpath('/labs/srslab/static_files/shared_apps/matlab_toolboxes/EEGLAB/eeglab2022.0/plugins/MicrostateAnalysis0.3')
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
study = 'DEFEND';
dir2load = '/labs/srslab/data_main/Microstates_vjp/binned/';
dir2save_csv = '/labs/srslab/data_main/Microstates_vjp/csvs/';
min_secs = 120;
epoch_length = 4;
sub_n = 25;
counter = 0;
plot_indiv_maps = 0;
% first create clusters maps for each person
for j = 1:sub_n
    subID = char(subIDs{j});
    EEG = pop_loadset( 'filename', [subID '_' study '_binned.set'],'filepath',dir2load);
    if size(EEG.ec,3) >= min_secs/epoch_length
        EEG.data = EEG.ec;
        EEG.trials = size(EEG.ec,3);
        EEG.epoch = EEG.epoch(1:EEG.trials);
        [EEG] = pop_FindMSTemplates(EEG, struct('MinClasses', 4, 'MaxClasses', 4, 'GFPPeaks', 1, 'IgnorePolarity', 1, 'MaxMaps', 1000, 'Restarts', 20, 'UseAAHC', 0), plot_indiv_maps, 0);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',[subID, ' indiv. microstates']);
        eeglab redraw
        counter = counter +1;
    end
end

indices = [1:counter];
%[EEG, com] = pop_CombMSTemplates(ALLEEG, indices, 0, 1, 'GrandMean');
[EEG, com] = pop_CombMSTemplates(ALLEEG, indices, 1, 1, 'DoMeans'); %doing 'DoMeans' instead of 'GrandMean' to better match da Cruz, but doesn't seem to make much of a difference
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, counter+1,'gui','off');
[ALLEEG com] = pop_SortMSTemplates(ALLEEG, indices, 0, 3);
grand_mean_set_idx = length(ALLEEG);
pop_QuantMSTemplates(ALLEEG, indices, 1, struct('b',5,'lambda',10,'PeakFit',1,'nClasses',4,'BControl',1), grand_mean_set_idx, [dir2save_csv,'microstate_results.csv']);



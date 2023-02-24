
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
study = 'DEFEND';
dir2load = '/labs/srslab/data_main/Microstates_vjp/Results/';
dir2save = '/labs/srslab/data_main/Microstates_vjp/binned/';
codes = ['11','12','13','21','22','23','50'];
bin_report=subIDs;
epoch_s = 4;
for j = 1:25 
    subID = char(subIDs{j});
    load([dir2load subID '_preprocess.mat']);
    EEG.ec = EEG.data(:,:,ismember(EEG.bin_vector,[21,22,23]));
    EEG.eo = EEG.data(:,:,ismember(EEG.bin_vector,[11,12,13]));
    %% then save out
    EEG = pop_saveset( EEG, 'filename',[subID '_' study '_binned.set'],'filepath', dir2save,'savemode','onefile');
    bin_report(j,2)={size(EEG.ec,3)*epoch_s};
    bin_report(j,3)={size(EEG.eo,3)*epoch_s};
end





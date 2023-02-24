% App beta version - Automatic Preprocessing Pipeline for EEG data
% This is the resting state alternative of App
% No guarantee that this works
% The current version is as of 16/10/2017
% By Janir Ramos da Cruz @ EPFL
% w/ modifications made by VJP to work with blocked data

% Directory management
%clc, clear all, close all
want2plot = 1;
codes = ['11','12','13','21','22','23','50'];
CurrDir = pwd;                                              % gets current directory
%eeglabpath = fileparts(which('eeglab.m'));                  % Getting eeglab path
addpath('/labs/srslab/data_main/Microstates_vjp/functions/')
addpath(genpath('APP_EEG/'))
% Paths to the data
allsubsDir2load = '/labs/srslab/data_main/Microstates_vjp/interp_recoded/';

% Info about your data
% Below are examples for the Schizophrenia Data
epoch = 4;                              % each epoch has 4 seconds for eeg analysis
n_ch = 64;                              % number of channels

% Info about the analysis
doICA = 1;                              % if wants to do ICA decomposition, set doICA = 1
d_Fs = 128;                             % desired sample rate in Hz
LowCutOffFreq = 1;                      % lower bound cut-off frequency
HighCutOffFreq = 40;                    % high bound cut-off frequency
filterOrder = [];                       % filterOrder - leave blank and determined by eeglab

% Creates the directory to save data for each group
SaveDir = '/labs/srslab/data_main/Microstates_vjp/Results/';
LoadDir = '/labs/srslab/data_main/Microstates_vjp/interp_recoded/';
study = 'DEFEND';

%% loop through subs
Subject_pool = subIDs;
for iSubject = 1:25
    % Initiate subject rejected data
    % Bad channels
    bad_channels = 0;
    % Bad trials
    bad_trials = 0;
    % Bad channels in each trial
    trial_bad_ch = 0;
    subID = char(Subject_pool{iSubject});
    % out load command here
    dataset = pop_loadset( 'filename', [subID '_' study '_interp_recoded.set'],'filepath',LoadDir);
    % Bandpass filter with 0 phase lag
    dataset = pop_eegfiltnew(dataset, LowCutOffFreq, HighCutOffFreq, filterOrder, 0, [], 0);

    % CleanLine to remove the 60 Hz and the harmonics
    dataset = pop_cleanline(dataset, 'Bandwidth',2,'ChanCompIndices',[1:dataset.nbchan],...
        'SignalType','Channels','ComputeSpectralPower',true,               ...
        'LineFrequencies',[60 120] ,'NormalizeSpectrum',false, ...
        'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,            ...
        'ScanForLines',true,'SmoothingFactor',100,'VerboseOutput',1,       ...
        'SlidingWinLength',4,'SlidingWinStep',4);

    % Downsample the data if needed
    if dataset.srate ~= d_Fs
        dataset = pop_resample(dataset,d_Fs);
    end

    % Calculate the EOG channels, increase the SNR
    %heog should be 133:134
    %veog should be 135:136
    veog = dataset.data(n_ch+2,:) - dataset.data(n_ch+1,:); % VEOG channel (left-right - looking at the subj)
    heog = dataset.data(n_ch+3,:) - dataset.data(n_ch+4,:); % HEOG channel (top - bottom)

    dataset.data(n_ch+1,:) = veog; % overwrite channels
    dataset.data(n_ch+2,:) = heog;
    clear veog heog % free some memory

    dataset = pop_select(dataset, 'channel', 1:n_ch+2); % remove channels 67 68
  
    %rename the EOG channels
    dataset.chanlocs(n_ch+1).labels = 'VEOG';
    dataset.chanlocs(n_ch+2).labels = 'HEOG';
    % extract EOG channels
    dataset_eog = pop_select(dataset, 'channel', n_ch+1:n_ch+2);
    % remove EOG channels from EEG data
    dataset = pop_select(dataset, 'channel', 1:n_ch);

    % Re-reference to the biweight mean of the channels
    % Because according to Biosemi the CMS electrode does not provide
    % the full 80 dB CMRR
    [avg_ch,~] = myBiweight(dataset.data');
    dataset.data = dataset.data - repmat(avg_ch,n_ch,1);
    clear avg_ch

    % Find bad channels
    % standard deviation
    [~,ch_std] = myBiweight(dataset.data);
    bad_ch_ind_1 = myFindOutliers(ch_std);
    % mean correlation of top 4 correlation coefficients
    ch_r_raw = corr(dataset.data');
    ch_r_raw = sort(ch_r_raw);
    ch_r = mean(ch_r_raw(end-4:end-1,:)); % top 4 correlation coefficients excluding self-correlation
    bad_ch_ind_2 = myFindOutliers(ch_r);
    % remove and interpolate bad channels
    bad_ch = unique([bad_ch_ind_1,bad_ch_ind_2]);
    if ~isempty(bad_ch)
        dataset = eeg_interp(dataset, bad_ch, 'spherical');
    end
    bad_channels = bad_channels + length(bad_ch); % the number of bad channels per subject
    %------------------------------------------------------------------------------------
    dataset.old_events = dataset.event;
    %epoch data in x seconds epochs
    dataset = eeg_regepochs( dataset, 'limits', [0 epoch], 'rmbase', NaN, 'recurrence', epoch);
    dataset = eeg_checkset( dataset );
    %if regeppochs gets rid of final event code then put it back in
    if ismember(50,unique([dataset.old_events.type])) & ~ismember('50',{dataset.event.type})
        [dataset.event(end).type] = '50';
        [dataset.epoch(end).eventtype] = '50';
    end

    dataset_eog = eeg_regepochs( dataset_eog, 'limits', [0 epoch], 'rmbase', NaN, 'recurrence', epoch);
    dataset_eog = eeg_checkset( dataset_eog );

    %give each epoch an event label
    eventtypes = {dataset.epoch.eventtype};
    block_start_idx = [];
    block_start_code = [];
    for jj = 1:length(eventtypes)
        single_event = eventtypes{jj};
        if contains(codes,single_event)
            single_code = single_event(~strcmp('X',single_event));
            if contains('0', single_event) && contains('11',single_event)
                single_code = '11';
            end
            block_start_idx = [block_start_idx,jj];
            block_start_code = [block_start_code, single_code];
        end
    end
    epochs_per_block = diff(block_start_idx);

    if ismember('50',block_start_code) || ismember('5',block_start_code)
        dataset.bin_vector = [str2double(repelem(block_start_code(1:end-1), epochs_per_block)),50];
    else
        epochs_per_block = [epochs_per_block, length(eventtypes)-block_start_idx(end)+1];
        dataset.bin_vector = str2double(repelem(block_start_code, epochs_per_block));
    end

    if size(dataset.data,3 ) - length(dataset.bin_vector) ==1
        dataset.bin_vector = [dataset.bin_vector,0];
    end 
    
    % Max amplitude difference
    amp_diffs = zeros(size(dataset.data,1),size(dataset.data,3));
    for iChan = 1:size(dataset.data,1)
        for itrial = 1:size(dataset.data,3)
            amp_diffs(iChan,itrial) = max(dataset.data(iChan,:,itrial)) - min(dataset.data(iChan,:,itrial));
        end
    end
    [epoch_amp_d,~] = myBiweight(amp_diffs');
    % Epoch variance or the mean GFP
    epoch_GFP = mean(squeeze(std(dataset.data,0,2)));
    % Epoch's mean deviation from channel means.
    [means,~] = myBiweight(dataset.data(:,:)); % channel mean for all epochs
    epoch_m_dev = zeros(1,size(dataset.data,3));
    for itrial = 1:size(dataset.data,3)
        epoch_m_dev(itrial) = mean(abs(squeeze(mean(dataset.data(:,:,itrial),2))' - means));
    end

    % Find the bad trials
    Rej_ep_amp_d = myFindOutliers(epoch_amp_d);
    Rej_ep_GFP = myFindOutliers(epoch_GFP);
    Rej_ep_mdev = myFindOutliers(epoch_m_dev);
    Rej_epoch = unique([Rej_ep_amp_d Rej_ep_GFP Rej_ep_mdev]);

    % Remove the bad trials and trim data so that there are no epochs after
    % the ending code ('50')
    extra_epochs = find(dataset.bin_vector==50):size(dataset.data,3);
    Rej_epoch  = [Rej_epoch,extra_epochs(2:end)];
    dataset.bin_vector = dataset.bin_vector(setdiff(1:length(dataset.bin_vector),Rej_epoch));
    dataset = pop_select(dataset,'notrial',Rej_epoch);
    dataset_eog = pop_select(dataset_eog,'notrial',Rej_epoch);
    % Count the bad trials
    bad_trials = bad_trials + length(Rej_epoch);

    if doICA == 1
        %------------------------------------------------------------------------------------
        % Perform ICA - SOBI
        %------------------------------------------------------------------------------------
        EEG = pop_runica(dataset, 'icatype', 'sobi', 'dataset',1, 'options',{});
        if isempty(EEG.icaact)
            disp('EEG.icaact not present. Recomputed from data.');
            if length(size(EEG.data))==3
                EEG.icaact = reshape(EEG.icaweights*EEG.icasphere*reshape(EEG.data,[size(EEG.data,1)...
                    size(EEG.data,2)*size(EEG.data,3)]),[size(EEG.data,1) size(EEG.data,2) size(EEG.data,3)]);
            else
                EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data;
            end
        end

        ncomp = length(EEG.icawinv); % number of components

        % Eye blinks and saccades detection by correlation with VEOG and HEOG
        VEOG = dataset_eog.data(1,:,:);
        VEOG = VEOG(:);
        HEOG = dataset_eog.data(2,:,:);
        HEOG = HEOG(:);
        ICs = EEG.icaact(:,:)';
        for ic = 1:size(ICs,2)
            corr_V(ic) = corr(ICs(:,ic),VEOG);
            corr_H(ic) = corr(ICs(:,ic),HEOG);
        end
        rej_V = myFindOutliers(corr_V);
        rej_H = myFindOutliers(corr_H);
        Rej_ic_eog = unique([rej_V,rej_H]);     % ICs containing blinks
        clear rej_V rej_H % free some memory

        % ICs with generics discontinutiy of spatial features
        topography = EEG.icawinv';  % topography of the IC weigths
        channel = EEG.chanlocs(EEG.icachansind');
        xpos=[channel.X];ypos=[channel.Y];zpos=[channel.Z];
        pos=[xpos',ypos',zpos'];
        gen_disc = zeros(1,size(ICs,2)); % generic discontinuity
        for ic = 1:ncomp
            aux = [];
            for el = 1:length(channel)-1
                P_el = pos(el,:); %position of current electrode
                d = pos - repmat(P_el,length(channel),1);
                dist = sqrt(sum((d.*d),2));
                [y,I] = sort(dist);
                rep_ch = I(2:11); % the 10 nearest channels to el
                weight_ch = exp(-y(2:11)); % respective weights, computed wrt distance
                aux = [aux abs(topography(ic,el)-mean(weight_ch.*topography(ic,rep_ch)'))];
                % difference between el and the average of 10 neighbor el
                % weighted according to weight
            end
            gen_disc(ic)=max(aux);
        end
        Rej_ic_gd = myFindOutliers(gen_disc); % ICs containing generic discontinuities
        clear channel aux pos topography % free some memory
        % Muscle activity usually has low autocorrelation of time course
        ncorrint =round(25/(1000/EEG.srate)); % number of samples for 25 ms lag
        for k = 1:ncomp
            y = EEG.icaact(k,:,:);
            yy = xcorr(mean(y,3),ncorrint,'coeff');
            autocorr(k) = yy(1);
        end
        Rej_muscle = myFindOutliers(autocorr);

        % Drop the ICs
        droplist = unique([Rej_ic_eog Rej_muscle Rej_ic_gd]); % save the number of the components dropped
        if ~isempty(droplist)
            EEG = pop_subcomp( EEG, droplist, 0); %drop components
        end
        % ------------------------------------------------------------------------------------
    else
        % when we do not use ICA
        EEG = dataset;
    end
    clear dataset
    % Check if still have some artifacts after IC -- focus on each
    % channel of each trial
    EEGtmp = EEG; % local copy
    for itrial = 1:EEG.trials
        EEGtmp.data = EEG.data(:,:,itrial);

        % Mean diff value
        [mean_diff,~]=myBiweight(diff(EEGtmp.data,[],2));
        % Variance of the channels
        [~,chan_std]=myBiweight(EEGtmp.data);

        % Find the outlier channels
        Rej_mean_diff = myFindOutliers(mean_diff);
        Rej_chan_std = myFindOutliers(chan_std);
        bad_ch_trial = intersect(Rej_mean_diff, Rej_chan_std); % 18/07/2017
        if ~isempty(bad_ch_trial)
            EEGtmp = eeg_interp(EEGtmp, bad_ch_trial, 'spherical');
        end
        % add the data to the dataset
        EEG.data(:,:,itrial) = EEGtmp.data;

        % count the number of bad channel in each trial
        trial_bad_ch = trial_bad_ch + length(bad_ch_trial);

    end
    clear EEGtmp % clear some memory

    %-------------------------------------------------------------------------------

    % Prepare the data to be stuck together
    % by using inverse hanning window
    xtmp = EEG.data(:,:);
    for ichan = 1:size(xtmp,1)
        % input arguments - data, sampling rate, fraction of the
        % data to apply the hanning, size of the data segment
        xtmp(ichan,:) = hann_intersection(xtmp(ichan,:), EEG.srate, 1/4, EEG.srate*epoch);
    end

    % copy data and check the set, the set check is already epoched
    % again but it's okay since already hanninged
    EEG.data = xtmp;

    % reference to the average
    EEG = pop_reref( EEG, []);
    
    if want2plot == 1
    figure; pop_spectopo(EEG, 1, [0      1992.1875], 'EEG' , 'percent', 100, 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
    ylims = ylim;
    colormap(redwhiteblue(ylims(1), ylims(2)));
    end

    save(strcat(SaveDir,subID,'_preprocess.mat'),'EEG',...
        'bad_trials','bad_channels','trial_bad_ch','droplist');

    if length(EEG.bin_vector) ~= size(EEG.data,3 )
        error('different lengths of epoch bins and data ... go figure out the problem')
    end

end

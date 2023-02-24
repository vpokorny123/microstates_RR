addpath('/labs/srslab/data_main/Microstates_vjp/functions')
SaveDir = '/labs/srslab/data_main/Microstates_vjp/interp_recoded/';
LoadDir = '/labs/srslab/data_staging/DEFEND_archive/';
chanlocs = readlocs('/labs/srslab/data_main/Microstates_vjp/pilot_scripts/matfiles/BioSemi_64_10-10.elp'); % reads the channels location with EEGLAB
codes = [11,21,12,22,13,23,50]; %these are the event codes denoting eyes open or eyes closed conditions
%% loop through subs
studies = {'DEFEND'};
for j = 1:length(studies)
    study = studies{j};
    for jj = 1:25
        subID = char(subIDs{jj});
        % Read in the BDF file using eeglab
        EEG = pop_biosig([LoadDir, '/', subID, '/', subID, '_RestingEEG.bdf'], 'ref',1,'importannot','off');
        if size(EEG.data,1) > 150 %check to see if data is BioSemi 128 if so ...
            %grab aux JIC we need them later
            EEG.origmontage = 'BioSemi_128';
            ears = pop_select(EEG,'channel',[129:130]);
            veog = pop_select(EEG,'channel',[131:132]);
            heog = pop_select(EEG,'channel',[133:134]);
            emg = pop_select(EEG,'channel',[135:138]);
            EEG = pop_select(EEG,'channel',[1:128]);
            %read in 128 montage info
            EEG = pop_chanedit(EEG,'load',{'/labs/srslab/static_files/shared_apps/matlab_toolboxes/ssk_eegtoolbox/ICAcleanEEG.v.1.3/montage/BioSemi_128_elecN.elp' 'filetype' 'autodetect'});
            %add ears back in
            EEG.data(end+1:end+2,:) = ears.data;
            EEG.chanlocs(end+1).label = 'L ear';
            EEG.chanlocs(end+1).label = 'R ear';
            % reref to ears
            EEG = pop_reref( EEG, [129 130]);
            %now interpolate
            EEG.data2interp = EEG.data(1:128,:);
            EEG.data = convert_128to64(EEG);
            EEG.nbchan = size(EEG.data,1);
            EEG.chanlocs = chanlocs;
        end
        %% lines 35- adjust event codes
        event_vec = [EEG.event.type];
        latencies = [EEG.event.latency];
        if ismember(768,event_vec)
            new_event_vec = event_vec-768;
            new_latencies = latencies;
            event_vec = new_event_vec;
        end
        if ismember(1023,event_vec)
            event_vec = event_vec-1023;
        end

        %the events with 1023/1024s are tricky, and require more attention
        if ismember(1024,event_vec)
            if any(event_vec==1024 & latencies<15000)
                starting_idx = event_vec==1024 & latencies<15000;
                starting_code = event_vec(starting_idx);
                starting_latency = latencies(starting_idx);
                if size(starting_code,1)>1
                    disp(['multiple valid starting points; choosing the last one but you should check; ' ...
                        'if starting point looks good press any key to continue'])
                    pause;
                    starting_code(size(starting_code,1),:)
                end
            elseif any(event_vec==0 & latencies<15000) %% if initial event doesn't exist then just use the initialization code
                starting_idx = event_vec==0 & latencies<15000;
                starting_code = event_vec(starting_idx);
                starting_latency = latencies(starting_idx);
            else
                error('could not find reasonable starting point go figure it out')
            end
            % now set block duration based on event latencies
            block_dur_s = 45;
            for jjj = 1:10
                if any(latencies>(EEG.srate*(block_dur_s+(jjj-1))) & latencies <EEG.srate*(block_dur_s+jjj))
                    block_dur = EEG.srate * (block_dur_s+(jjj-1));
                end
            end
            % so now we rebuild our event codes a priori using our
            % starting point and block durations
            new_latencies = zeros(1,7);
            for jjj = 1:length(new_latencies)
                new_latencies(jjj) = starting_latency + (block_dur*(jjj-1));
            end
            new_event_vec = codes;
            %checking to see if data collector terminated recording before task
            %ended
            min_run_time_mins = 4.5;
            min_run_time = min_run_time_mins*60*EEG.srate;
            if max(latencies)< min_run_time && size(EEG.data,2)>min_run_time
                disp("looks like the recording is shorter than it should be" + ...
                    "you should check to see why this might be. If not a problem" + ...
                    "press any key to continue and ending event code will be added" + ...
                    "to final timepoint of recording")
                latencies(end+1) = size(EEG.data,2);
            end
            % now compare old and new latencies and copy over the "true
            % latency"  that is nearest to the estimated latency
            for jjj = 2:length(new_latencies)
                single_new_latency = new_latencies(jjj);
                latency_diff = single_new_latency-latencies;
                [~,closest_latency_idx] = min(abs(latency_diff));
                new_latencies(jjj) = latencies(closest_latency_idx);
            end
        end
        cell_latencies = num2cell(new_latencies');
        cell_events= num2cell(new_event_vec');
        [EEG.event(length(cell_events)+1:end)] = [];
        [EEG.event.type] = cell_events{:};
        [EEG.event.latency] = cell_latencies{:};
        %finally add veog and heog channels back in
        EEG.data = [EEG.data; veog.data; heog.data];
        EEG.chanlocs(end+1).labels = 'Hi VEOG'; %hi 65
        EEG.chanlocs(end+1).labels = 'Lo VEOG'; %lo 66
        EEG.chanlocs(end+1).labels = 'R HEOG'; % R 67
        EEG.chanlocs(end+1).labels = 'L HEOG'; % L 68
        %get rid of data2interp to save space
        EEG = rmfield(EEG,'data2interp');
        % then save out
        EEG = pop_saveset( EEG, 'filename',[subID '_' study '_interp_recoded.set'],'filepath', SaveDir,'savemode','onefile');
    end
end

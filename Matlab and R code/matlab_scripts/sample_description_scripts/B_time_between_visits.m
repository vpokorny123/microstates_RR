%directories of interest
load('/labs/srslab/data_main/Microstates_vjp/sample_description_scripts/matfiles/subject_counts.mat')
direcs = {'/labs/srslab/data_staging/NPG_archive/',
    '/labs/srslab/data_staging/GTF_archive/',
    '/labs/srslab/data_staging/PENS_EEG/'};
studies = {'NPG','GTF','PENS'}; % make sure this list is in the same order as the directory list above
end_names = {'.bdf','.bdf','.vhdr'};


subIDs = at_least_twos;
day_report = subIDs;
for j = 1:length(subIDs)
    tic
    subID = char(subIDs{j});
    dates = [];
    for jj = 1:length(studies)
        study = char(studies{jj});
        %search through direcs
        study_id_list= resting_list(strcmp(resting_list(:,1),subID),2);
        if any(strcmp(study_id_list,study))
            end_name = char(end_names{jj});
            direc = [char(direcs{jj}),subID,'/'];
            file_list = dir([direc,'*',end_name]);
            unzipped_file_list = {file_list(:).name};
            resting_idx = find(contains(unzipped_file_list,'esting'));
            resting_file = char(unzipped_file_list(resting_idx));
            if size(resting_file,1)>1
                resting_file = resting_file(1,:);
            end
            if strcmp(end_name, '.bdf')
                EEG = pop_biosig([direc,resting_file]);
                dates=[dates;EEG.etc.T0];
            elseif strcmp(end_name,'.vhdr')
                %EEGLABS brainvision loader doesn't give you dates automatically...
                fileID = fopen([direc,subID '_RestingEEG.vmrk']);
                bv_info = textscan(fileID,'%s','Delimiter',',');
                raw_date = bv_info{1}{20};
                clean_date = [str2double(raw_date(1:4)),str2double(raw_date(5:6)), ...
                               str2double(raw_date(7:8)),str2double(raw_date(9:10)), ...
                               str2double(raw_date(11:12)),str2double(raw_date(13:14))];
                dates = [dates;clean_date];
            end  
        end
    end
    if size(dates,2) == 3
        start_date = dates(1,1:3);
        end_date = dates(3,1:3);
    else
        start_date = dates(1,1:3);
        end_date = dates(2,1:3);
    end

    day_report(j,2) = {daysact(datetime(start_date),datetime(end_date))};
    day_report(j,3) = {yearfrac(datetime(start_date),datetime(end_date))};
    disp(['done with subj #',num2str(j)])
    toc
end

 mean_in_days = mean([day_report{:,2}]);
 mean_in_years =  mean([day_report{:,3}]);
 mean_in_years
 std([day_report{:,3}])
 histogram([day_report{:,3}])

 
 
      
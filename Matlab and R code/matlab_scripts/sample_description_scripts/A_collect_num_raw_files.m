%directories of interest
direcs = {'/labs/srslab/data_staging/GTF_archive/',
    '/labs/srslab/data_staging/NPG_archive/',
    '/labs/srslab/data_staging/PENS_EEG/'};
studies = {'GTF','NPG','PENS'}; % make sure this list is in the same order as the directory list above

all_storage = {};
resting_storage = {};
resting_studyid_storage = [];
all_studyid_storage = [];
for j = 1:length(direcs)
    directory = direcs{j};
    counter = dir([directory, '*00' ]);
    all_storage = [all_storage;{counter.name}'];
    all_studyid_storage = [all_studyid_storage; cellstr(repmat(studies{j},length({counter.name}),1))];
    for jj = 1:length({counter.name})
        sub_files = dir([directory,counter(jj).name]);
        if any(contains({sub_files.name},'esting'))
            resting_storage = [resting_storage; counter(jj).name];
            if contains(directory,'GTF')
                resting_studyid_storage = [resting_studyid_storage;cellstr('GTF')];
            elseif contains(directory,'NPG')
                resting_studyid_storage = [resting_studyid_storage;cellstr('NPG')];
            elseif contains(directory,'PENS')
                resting_studyid_storage = [resting_studyid_storage;cellstr('PENS')];
            end
        end  
    end
    
end



resting_list = [resting_storage,resting_studyid_storage];
all_list = [all_storage,all_studyid_storage];
 
 %% look into folks without resting files to make sure they don't have just weird names
 for j=1:length(studies)
 missing_rest.(studies{j}) = setdiff(all_list(strcmp(all_list(:,2),studies{j})),resting_list(strcmp(resting_list(:,2),studies{j})));
 end
 
 %% look into how many repeats we have
 [gc,grps] = groupcounts(resting_list(:,1));
 singles = grps(gc ==1);
 dubs = grps(gc ==2);
 trips = grps(gc ==3);
 at_least_ones = length(singles)+length(dubs)+length(trips);
 at_least_ones_alt = grps(gc>=1);
 at_least_twos = grps(gc>=2);
 save('/labs/srslab/data_main/Microstates_vjp/sample_description_scripts/matfiles/subject_counts.mat','at_least_ones','at_least_twos','trips','dubs','singles','resting_list')
 
 writetable(cell2table(resting_list),'/labs/srslab/data_main/Microstates_vjp/csvs/resting_subs_raw.csv');
 

 
 
      
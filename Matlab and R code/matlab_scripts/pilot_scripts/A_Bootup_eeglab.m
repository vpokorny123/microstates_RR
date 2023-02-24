% preprocess resting files from DEFEND sutdy of mTBI and PTSD 
% for building and testing processing pipeline
addpath('/labs/srslab/static_files/shared_apps/matlab_toolboxes/EEGLAB/eeglab2021.1/')
eeglab

% Subjects BDF files
Subject_dir = '/labs/srslab/data_staging/DEFEND_archive/';
Subject_BDF = dir(Subject_dir);
Subject_pool = regexpi({Subject_BDF.name},'[0-9]{4}','match');
Subject_pool = Subject_pool(~cellfun('isempty',Subject_pool))';
subIDs = Subject_pool;

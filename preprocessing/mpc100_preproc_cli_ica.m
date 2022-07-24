function mpc100_preproc_cli_ica(sub_n)
% sub_n = 32; % for test

%% Basic Setting

addpath(genpath('/home/donghee/Documents/resources/github/canlab')) % CANLAB CORE
addpath(genpath('/home/donghee/Documents/resources/github/cocoanlab')) % COCOAN CORE
addpath(genpath('/home/donghee/Documents/resources/spm12')) % SPM
addpath('/cocoanlab2/GPU3_sync/projects/MPC/MPC_100/sync/scripts/imaging')


%% SETUP(1): Set run number, task names, and number of disdaq (Usually same for the every subject)

func_run_nums = 1:10;
run_n = length(func_run_nums);
func_tasks = {'resting', 'nomovie', 'movie', 'movie', 'movie', 'movie', 'movie', 'movie', 'nomovie', 'caps'};
disdaq_n = repmat(18, 1, run_n); %(number of TR, 1, number of run);

%% SETUP(2): Making directory and subject code

basedir = '/cocoanlab2/GPU3_sync/data/MPC/MPC_100';

ica_aroma_basedir = '/home/donghee/ICA-AROMA-MPC';

study_imaging_dir = fullfile(basedir, 'imaging');


subj_idx = sub_n;

projName = 'mpc'; % project name
[preproc_subject_dir, subject_code] = mpc100_make_subject_dir_code(study_imaging_dir, projName, subj_idx);

num_sub = length(subject_code);

% once compeleting one subject, then do it for all other subjects in HPC
% preproc_subject_dir = preproc_subject_dir{2};
% subject_code = subject_code{1,2};
% num_sub = length(subject_code);


%% A. BIDS dicom to nifti =================================================
% PART A: --------------------
% 1. dicom to nifti: bids
% 2. bids validation
% ----------------------------

% You can make these as a loop for multiple subjects

%% A-1. Make directories ('raw' directory will be made)
for i=1:num_sub
    humanfmri_a1_make_directories(subject_code{1,i}, study_imaging_dir, func_run_nums, func_tasks);
end

% After this command, you have to move the directories that contain dicom files 
% into the corresponding directories. 

%% Move Directories ('dicom_from_scanner' -> 'raw)
addpath('/media/das/cocoanlab Dropbox/projects/MPC/MPC_100/sync/scripts/imaging')

for i = 1:num_sub
    mpc100_move_dicom_to_raw(subject_code{1,i}, study_imaging_dir, run_n, 'copy');
end


%% remove the subjects who were scanned two times or had some issues in the scan time.
% MPC002(idx 1) MPC009(idx 6) MPC016(idx 13) MPC030(idx 24) MPC038(idx 30) MPC052(idx 43)
% MPC054(idx 44) MPC055(idx 45) MPC066(idx 52) MPC067(idx 53) MPC078(idx 62) 
% MPC090(idx 74) MPC098(idx 79) MPC107(idx 88)


 subj_idx_rm = [2:5, 7:12, 14:23, 25:29, 31:42, 46:51, 54:61, 63:73, 75:78, 80:87, 89:91];
 
subject_code_rm = subject_code(subj_idx_rm);
 
num_sub_rm = length(subject_code_rm);

%% A-2. Dicom to nifti: structural and functional 
d=datetime('now');
 
for i=67:num_sub_rm
    % A-2. Dicom to nifti: anat(T1)
    
    humanfmri_a2_structural_dicom2nifti_bids(subject_code_rm{1,i}, study_imaging_dir);
    
    % A-3. Dicom to nifti: functional(Run 1~10)
    
    humanfmri_a3_functional_dicom2nifti_bids(subject_code_rm{1,i}, study_imaging_dir, disdaq_n, 'no_check_disdaq');
    
    % A-4. Dicom to nifti: fmap(Distortion correction)
    
    humanfmri_a4_fieldmap_dicom2nifti_bids(subject_code_rm{1,i}, study_imaging_dir);
    
    d=[d datetime('now')];
    
end

%% for one subject

% A-2. Dicom to nifti: anat(T1)
%humanfmri_a2_structural_dicom2nifti_bids(subject_code, study_imaging_dir);

% A-3. Dicom to nifti: functional(Run 1~10)

%humanfmri_a3_functional_dicom2nifti_bids(subject_code, study_imaging_dir, disdaq_n);
%humanfmri_a4_fieldmap_dicom2nifti_bids(subject_code, study_imaging_dir);


%% Done with A Part =======================================================

% You can use the following tools

% BIDS-validator: http://incf.github.io/bids-validator/ in chrome
% Fmriprep (preprocessing tool): see http://fmriprep.readthedocs.io/en/stable/

%% B. COCOANLAB PREPROC

% PART B: --------------------
% 3. disdaq & visualization/qc (canlab)
% 4. motion correction (realignment) 
% 5. EPI normalization 
% 6. Smoothing
% 7. ICA-AROMA 
% ----------------------------

%% PART 3 
d=datetime('now')

% B-1. Preproc directories
humanfmri_b1_preproc_directories(subject_code, study_imaging_dir); %'forced_save', 'no_save'

%% B-2. Implicit mask and save means

% preproc_subject_dir 
humanfmri_b2_functional_implicitmask_savemean(preproc_subject_dir);

%% B-3. Spike id
humanfmri_b3_spike_id(preproc_subject_dir);

%% B-4. Slice timing correction if needed: You can skip this if TR is short enough

% tr = .46;
% mbf = 8;
% 
% humanfmri_b4_slice_timing(preproc_subject_dir, tr, mbf);

%% B-5. Motion correction

use_st_corrected_data = false;
use_sbref = true;
humanfmri_b5_motion_correction(preproc_subject_dir, use_st_corrected_data, use_sbref);

%% B-6. distortion correction    
epi_enc_dir = 'ap';
use_sbref = true;
humanfmri_b6_distortion_correction(preproc_subject_dir, epi_enc_dir, use_sbref, 'run_num', func_run_nums)

%% PART 4
% B-7. coregistration (spm_check_registration.m)

use_sbref = true;
%humanfmri_b7_coregistration(preproc_subject_dir, use_sbref);
humanfmri_b7_coregistration(preproc_subject_dir, use_sbref, 'no_check_reg');

%% B-8-1. T1 Normalization

use_sbref = true;
%humanfmri_b8_normalization(preproc_subject_dir, use_sbref);
humanfmri_b8_normalization(preproc_subject_dir, use_sbref, 'no_check_reg');

%% B-9. Smoothing

humanfmri_b9_smoothing(preproc_subject_dir);

%% B-10. ICA-AROMA

n_dim = 200;
%n_thread = 1;
humanfmri_b10_ICA_AROMA(preproc_subject_dir, 'ica_aroma_dir', ica_aroma_basedir, 'dim', n_dim)

fprintf('done. \n')

%% Step C: Check Framewise Displacement and make Nuisance Regressors

% humanfmri_c1_move_clean_files
%% C-2
humanfmri_c2_get_framewise_displacement(preproc_subject_dir)
% 
%% C-3
humanfmri_c3_make_nuisance_regressors(preproc_subject_dir)
% %make_nuisance_regressors(PREPROC,'regressors',{'24Move','Spike','WM_CSF'})
% %make_nuisance_regressors(PREPROC,'img','swr_func_bold_files')
%% DONE %%
fprintf('done. \n')
datetime('now')
end

% MPC temporal SNR (tSNR)
% modified version of fast_201204_tsnr.m

clc;
clear;

%% Path set up
addpath('/media/das/cocoanlab Dropbox/projects/MPC/MPC_100/sync/function')

wh_location = 'gpu3';

[basedir, datdir, projdir, behavdir] = set_path_MPC(wh_location);
sync_dir =  fullfile(projdir, 'sync');
img_dir = fullfile(datdir, 'imaging', 'preprocessed');

%% Load Dataset D

dataset_file = filenames(fullfile(behavdir, 'MPC_*_dataset_behavioral*.mat'));
dataset_file = dataset_file{numel(dataset_file)};

if isempty(dataset_file)
    error('There is no Dataset file')
end

load(char(dataset_file)); % variable name : D


%% subject list
% sub004 has no 'fmri_start_time
sub_list = [2 3 4 5 8 9 10 11 12 13 ...
            14 15 16 17 18 19 21 22 24 25 ...
            26 28 29 30 31 32 34 35 37 38 ... 
            39 41 42 43 44 45 46 47 48 49 ...
            50 51 52 54 55 59 60 61 63 64 ...
            65 66 67 68 69 70 71 73 74 75 ...
            76 78 79 80 81 82 83 84 85 86 ...
            87 88 89 90 91 92 93 97 98 99 ...
            100 101 102 103 104 105 106 107 108 110 ...
            111 112 114 116 118 119 120 121 123 124 ...
            125 126 127 128 130 131 132 133 134 135 ...
            136 137 138 139 141 144 145 146 147 150 ...
            151 152 153 154];

sub_dir = {};

for sub_i = 1:numel(sub_list)
    sub_dir{sub_i} = sprintf('sub-mpc%03d', sub_list(sub_i)); 
end

subject_dir = fullfile(img_dir, sub_dir);


%% Run-level tSNR map ( in the case without trial-level tSNR map )

run_nums = [2:9]; % Heat run only

for sub_i = 1:numel(sub_list)
    tic;
    
    % directory
    
    tsnr_dir = fullfile(subject_dir{sub_i},'tsnr');
    if ~exist(tsnr_dir, 'dir'), mkdir(tsnr_dir); end
    
    tsnr_run_level_dir = fullfile(tsnr_dir, 'run-level'); % run-level directory
    if ~exist(tsnr_run_level_dir, 'dir'), mkdir(tsnr_run_level_dir); end
    
    tsnr_subj_level_dir = fullfile(tsnr_dir, 'subject-level'); % subject-level directory
    if ~exist(tsnr_subj_level_dir, 'dir'), mkdir(tsnr_subj_level_dir); end
        
    func_nomovie_dir = filenames(fullfile(subject_dir{sub_i}, 'func','r*-mpc*_task-nomovie_run-*_bold.nii'));
    func_movie_dir = filenames(fullfile(subject_dir{sub_i}, 'func','r*-mpc*_task-movie_run-*_bold.nii'));
    func_dir = vertcat(func_nomovie_dir, func_movie_dir);
    
    load(fullfile(subject_dir{sub_i}, 'PREPROC.mat'));
    [pathstr, fname, ext] = fileparts(PREPROC.coreg_anat_file);
    
    anat_dir = fullfile(subject_dir{sub_i}, 'anat');
    
    deform_img = fullfile(anat_dir, ['y_' fname '.nii']);
    
    
    % implicit mask
    
    humanfmri_b2_functional_implicitmask_savemean_byeol(subject_dir(sub_i));
    source = fullfile(subject_dir{sub_i}, 'implicit_mask_99.nii');
    mask_img = fullfile(tsnr_dir, 'implicit_mask_99.nii');
    movefile(source, mask_img);
    
    tic;
    for run_i = 1:numel(run_nums)
        %target_img = fullfile(tsnr_dir, ['tsnr_run0' num2str(run_nums(run_i)) '.nii']);
        target_img = fullfile(tsnr_run_level_dir, sprintf('tsnr_run%02d.nii', run_nums(run_i)));
       
        % tsnr
        dat = fmri_data(func_dir{run_i}, mask_img);
        mean_dat = mean(dat);
        tsnr = dat;
        tsnr.dat = mean_dat.dat./std(dat.dat,1,2);
        tsnr.fullpath = target_img;
        write(tsnr);
        
        % normalization
        clear matlabbatch
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deform_img};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {target_img};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78  -112   -70
            78    76    85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_';
        
        spm('defaults','fmri');
        spm_jobman('initcfg');
        spm_jobman('run', {matlabbatch});
    end
    toc;
end

%% Subject-level tSNR map
for sub_i = 1:numel(sub_list)
    tsnr_run_level_dir = fullfile(subject_dir{sub_i},'tsnr','run-level');
    tsnr_subj_level_dir=fullfile(subject_dir{sub_i},'tsnr','subject-level');
    
    target_img = filenames(fullfile(tsnr_run_level_dir, 'w_tsnr_run*.nii'));
    if numel(target_img) == 8
        mdat = mean(fmri_data(target_img));
        mdat.fullpath = fullfile(tsnr_subj_level_dir, sprintf('mean_tsnr_%03d.nii', sub_list(sub_i)));
        write(mdat)
    else
        break
    end
end

%% Group-level tSNR map: (1) mean map (2)t-value map
subject_dir = fullfile(img_dir, sub_dir);

gray_matter_mask = which('gray_matter_mask.img');

for sub_i = 1:numel(sub_list)
    
    tsnr_dir = fullfile(subject_dir{sub_i},'tsnr');
    data = fmri_data(filenames(fullfile(tsnr_dir, 'subject-level','mean_tsnr_*.nii')), gray_matter_mask);
    if sub_i == 1
        target_img = data;
    else
        target_img.dat(:,sub_i) = data.dat;
    end
end

map_dir = fullfile(sync_dir, 'results', 'maps');
if ~exist(map_dir, 'dir'), mkdir(map_dir); end

% mean map
mdat = mean(target_img);
mdat.fullpath = fullfile(map_dir, 'mpc100_heat_group-level_tsnr_mean.nii');
write(mdat);

%t-value map
tdat = ttest(target_img, .05, 'bfr');
tdat.volInfo = data.volInfo;
tdat.fullpath = fullfile(map_dir, 'mpc100_heat_group-level_tsnr_tvalue.nii');
write(tdat);

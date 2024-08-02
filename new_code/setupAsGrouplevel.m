% This will load in all the voxelwise connectivity matrices and store them
% from all patients in 1 variable. 

clear; 

%% Setup the data for predictions - 4 voxel cubes
load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/HCP_IQ.mat');

cpath = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/voxelwise_connmats/';
cfiles = dir([cpath 'HCP_*_4vox.mat']);

ix=0;
for f = 1:length(cfiles)

    load([cfiles(f).folder '/' cfiles(f).name]);
    As(:,:,f) = connmat;
    clear connmat

    % pull IQ data from only the participants with connectomes
    for ff = 1:length(HCP_IQ)
        samesub = strcmp(cfiles(f).name(5:10),num2str(HCP_IQ(ff,1)));
        
        if samesub == 1
            ix=ix+1;
            Y(ix,1) = HCP_IQ(ff,2);
            fprintf('%d ',ix);
        end
    end
end
fprintf('\n');
outf = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/HCP_4voxel_connmats_all_subs.mat';
save(outf,"As",'-v7.3');

%% Do it again for the 3 voxel cubes
clear;

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/HCP_IQ.mat');

cpath = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/voxelwise_connmats/';
cfiles = dir([cpath 'HCP_*_3vox.mat']);

ix=0;
for f = 1:length(cfiles)

    load([cfiles(f).folder '/' cfiles(f).name]);

    % pull IQ data from only the participants with connectomes
    for ff = 1:length(HCP_IQ)
        samesub = strcmp(cfiles(f).name(5:10),num2str(HCP_IQ(ff,1)));
        
        if samesub == 1
            ix=ix+1;

            As(:,:,ix) = connmat;
            clear connmat

            Y(ix,1) = HCP_IQ(ff,2);
            fprintf('%d ',ix);

            % save every 100 subs in case it crashes
            if rem(ix/100,1) == 0
                outf = sprintf('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/HCP_3voxel_connmats_%d_to_%d_subs.mat',f-99,f);
                save(outf,"As",'-v7.3');
                clear outf As
                ix=0;
            end
        end
    end
end


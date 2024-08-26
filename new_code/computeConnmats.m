% This script computes the voxelwise connectivity matrices on just one
% hemisphere at a time (no cerebellum or subcortical)

clc; clear; close all;

addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

data_path = '/mnt/dustin/data/REST_LR/preprocessd_output/';
save_path = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Voxel_connmats/';

all_d = dir([data_path, '*_REST_LR_GSR_preprocessed_output.nii.gz']);
subvols_path = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Subvolumes_info/';

no_sub = length(all_d);

ix=0;
side = {'R','L'};

for f = 1:length(all_d)

    isvox = strfind(all_d(f).name,'REST_LR_GSR_preprocessed_output.nii.gz');

    if ~isempty(isvox)
        ix=ix+1;
        file = sprintf('%s%s',data_path,all_d(f).name);
        data = load_untouch_nii(file);

        img = data.img;

        % Loop through both sides
        for sd = 1:length(side)
            subvols_path_temp = sprintf('%sSubvolumes_%s/',subvols_path,side{sd});
            subvols_temp = dir([subvols_path_temp '*.mat']);

            % loop through each subvolume
            for sv = 1:length(subvols_temp)

                load([subvols_temp(sv).folder '/' subvols_temp(sv).name]);

                img_sub = zeros(size(img));

                for v = 1:length(subvol)
                    img_sub(subvol(v,1),subvol(v,2),subvol(v,3),:) = img(subvol(v,1),subvol(v,2),subvol(v,3),:);
                end

                dvec = reshape(img_sub,[(size(img_sub,1)*size(img_sub,2)*size(img_sub,3)),size(img_sub,4)]);

                %remove rows of zeros
                dvec(~any(dvec,2),:) = [];
                dvec = double(dvec);

                connmat(:,:) = double(corr(dvec(:,:)',dvec(:,:)'));
                connmat = normalize(connmat);

                fnum = subvols_temp(sv).name(15:16);

                outpath = sprintf('%s/Subvol_%s_%s',save_path,fnum,side{sd});  % ,side{s});

                if ~exist("outpath",'dir')
                    mkdir(outpath);
                end

                outfile = [outpath '/' all_d(f).name(1:11)];
                save(outfile,'connmat','-v7.3');

                clear dvec img_sub connmat outfile outpath fnum
            end
            clear subvols_path_temp subvols_temp
        end
    end
end

fprintf('\n');

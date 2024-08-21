% This script computes the voxelwise connectivity matrices on just one
% hemisphere at a time (no cerebellum or subcortical)

clc; clear; close all;

addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

data_path = '/mnt/dustin/data/REST_LR/preprocessd_output/';
save_path = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/';

all_d = dir([data_path, '*_REST_LR_GSR_preprocessed_output.nii.gz']);

no_sub = length(all_d);

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/JointParcelPredict-main/new_code/HCP_overlap_mask.mat');


%% 
ix=0;
side = {'Right','Left'};

for f = 1:length(all_d)

    isvox = strfind(all_d(f).name,'REST_LR_GSR_preprocessed_output.nii.gz');

    if ~isempty(isvox)
        ix=ix+1;
        file = sprintf('%s%s',data_path,all_d(f).name);
        data = load_untouch_nii(file);

        img = data.img;

        % populate this patient's image with only voxels in the masks
        img_cleaned = zeros(size(img));

        for v = 1:size(x_pos)
            img_cleaned(x_pos(v),y_pos(v),z_pos(v),:) = img(x_pos(v),y_pos(v),z_pos(v),:);
        end

        dvec = reshape(img_cleaned,[(size(img_cleaned,1)*size(img_cleaned,2)*size(img_cleaned,3)),size(img_cleaned,4)]);

        %remove rows of zeros
        dvec(~any(dvec,2),:) = [];
        dvec = double(dvec);

        connmat(:,:) = double(corr(dvec(:,:)',dvec(:,:)'));
        connmat = normalize(connmat);

        outfile = sprintf('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/compare_to_shen368/HCP_%s.mat',...
            all_d(ix).name(1:11));  % ,side{s});

        save(outfile,'connmat','-v7.3');
        clear dvec img_cleaned connmat outfile % x_pos y_pos z_pos
%         end
        fprintf('%d ',ix);
    end
end

fprintf('\n');

% This script was used to create the HCP overlap mask

clc; clear; close all;

addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

data_path = '/mnt/dustin/data/REST_LR/preprocessd_output/';
save_path = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/';

all_d = dir([data_path, '*_REST_LR_GSR_preprocessed_output.nii.gz']);

no_sub = length(all_d);

for i = 1: no_sub
    nii = load_untouch_nii([all_d(i).folder, '/', all_d(i).name]); 
    img = nii.img;
    img_norm = sum(img.^2, 4);
    
    if( i==1)
        dim = size(img_norm);
        mask_all = zeros([dim, no_sub]);
    end
    mask_all(:,:,:,i)= double(img_norm>0);
end

mask_all= uint8(mask_all);

save([save_path, 'HCP_overlap_mask'], 'mask_all');

non_zero_coords = true(size(mask_all,1:3));

for i = 1:size(mask_all, 1)
    for j = 1:size(mask_all, 2)
        for k = 1:size(mask_all, 3)
            % Check if the pixel is non-zero across all images
            if any(mask_all(i,j,k,:) == 0)
                non_zero_coords(i,j,k) = false;
            end
        end
    end
end

% Index the non-zero coordinates
[x_pos, y_pos, z_pos] = ind2sub(size(mask_all,1:3), find(non_zero_coords));

clearvars -excep non_zero_coords x_pos y_pos z_pos mask_all save_path

save([save_path '/HCP_overlap_mask.mat']);

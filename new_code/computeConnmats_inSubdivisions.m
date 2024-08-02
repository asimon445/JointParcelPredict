% This script will prepare the voxel time series data to be inputted into
% the joint parcellation-prediction algorithm. 
%
% It will do so by systematically dividing up the brain into smaller cubes
% (14x14x14 voxels) and compute the voxelwise connectivity matrix for each
% cube. This needs to be done because the full voxelwise connectivity
% matrices are too big for the server to handle. 
% 
% Each division will be done hemisphere at a time, and will not include
% cerebellar or subcortical voxels.
%
% The output will contain a 2744 x 2744 x n (n=number of subjects)
% voxelwise connectivity matrix for each subregion (2744 is the number of 

clc; clear; close all;

addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

data_path = '/mnt/dustin/data/REST_LR/preprocessd_output/';
save_path = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/';

all_d = dir([data_path, '*_REST_LR_GSR_preprocessed_output.nii.gz']);

no_sub = length(all_d);

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/JointParcelPredict-main/new_code/HCP_overlap_mask.mat');
load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/JointParcelPredict-main/new_code/Shen368_10network.mat');

Shen368 = load_untouch_nii('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/JointParcelPredict-main/new_code/Shen_368_2mm.nii.gz');
Shen368 = Shen368.img;

coords_R = logical(false(91,109,91));
coords_L = logical(false(91,109,91));

nodeIx = find(Shen368_10network(:,2) >= 1 & Shen368_10network(:,2) <= 8);

for i = 1:size(Shen368, 1)
    for j = 1:size(Shen368, 2)
        for k = 1:size(Shen368, 3)
            % make sure that this voxel is in the HCP overlap mask
            if non_zero_coords(i,j,k) == 1

                thisvox_node = Shen368(i,j,k);

                % check if this voxel is in any of the 8 networks we care about
                if ismember(thisvox_node,nodeIx)
                    % check if this voxel is in the right or left hemisphere
                    sideIx = Shen368_10network(thisvox_node,3);
                    
                    if sideIx == 1
                        coords_R(i,j,k) = true;
                    elseif sideIx == 2
                        coords_L(i,j,k) = true;
                    end
                    clear sideIx
                end
            end
            clear thisvox_node
        end
    end
end

[x_pos_R, y_pos_R, z_pos_R] = ind2sub(size(mask_all,1:3), find(coords_R));
[x_pos_L, y_pos_L, z_pos_L] = ind2sub(size(mask_all,1:3), find(coords_L));

clear side_ix coords_reduced nodeIx x_pos y_pos z_pos


num_non_zero = length(x_pos);
neighbormat = zeros(num_non_zero,num_non_zero);

for i = 1:num_non_zero
    for j = 1:num_non_zero
        if i ~= j
            % Check if positions are neighbors (26 possibilities)
            if abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 0 && abs(z_pos(i) - z_pos(j)) == 0
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 0 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 0
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 0 && abs(y_pos(i) - y_pos(j)) == 0 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 0
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 0 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 0 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            end
        end
    end
end

total_neighbors = sum(neighbormat, 2);

save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/compare_to_shen368/neighbormat.mat','neighbormat');
save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/compare_to_shen368/total_neighbors.mat','total_neighbors');



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

        % for s = 1:2

        % populate this patient's image with only voxels in the masks
        img_cleaned = zeros(size(img));
% 
%         if s == 1
%             x_pos = x_pos_R;
%             y_pos = y_pos_R;
%             z_pos = z_pos_R;
%         else
%             x_pos = x_pos_L;
%             y_pos = y_pos_L;
%             z_pos = z_pos_L;
%         end

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
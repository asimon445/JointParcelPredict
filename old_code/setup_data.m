% This script will load in the voxelwise time series data from the first
% resting state scan in the transdiagnostic dataset and set it up to be run
% through the joint parcellation-prediction algorithm.
%
% For development purposes, the size of the data needs to be reduced. So,
% instead of computing ~150000 x 150000 connectivity matrices for each
% individual, the data will be downsampled to 1000x1000

clc; clear; close all;

files = dir('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Data/');
addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/sublist_and_bsiDep.mat');

% nvoxels = 1000;   % This is the number of voxels we'll downsample to

% load in fMRI data and compute voxel-wise connectivity matrices
idx = 0;
for f = 1:length(files)
    isnii = strfind(files(f).name,'nii.gz');
    if ~isempty(isnii)
    
        for ff = 1:length(newsublist)

            submatch = strfind(lower(files(f).name),lower(newsublist{ff,1}));

            if ~isempty(submatch)
                idx=idx+1;
                final_sublist{idx,1} = newsublist{ff,1};
                final_BSI_dep(idx,1) = new_bsi(ff,1);

                data = load_untouch_nii([files(f).folder '/' files(f).name]);
                Y = data.img;
                Y(~any(Y,1:3),:) = [];

                % Define the size of the cubes
                cube_size = [3, 3, 3, 1];   % 1 in position 4 is just what the number of samples will be divided by in the time domain
                
                % Calculate the size of the output matrix
                output_size = round(size(Y) ./ cube_size);

                % Initialize the output matrix
                smallerY = zeros(output_size,'int32');

                % Loop through each cube
                for i = 1:output_size(1)
                    for j = 1:output_size(2)
                        for k = 1:output_size(3)
                            % Calculate the indices of the current cube
                            cube_indices_x = (i-1)*cube_size(1) + (1:cube_size(1));
                            cube_indices_y = (j-1)*cube_size(2) + (1:cube_size(2));
                            cube_indices_z = (k-1)*cube_size(3) + (1:cube_size(3));

                            % Extract the cube from the input matrix
                            cube = Y(cube_indices_x, cube_indices_y, cube_indices_z, :);

                            % Average the cube along the first three dimensions
                            averaged_cube = mean(cube, [1, 2, 3]);

                            % Assign the averaged cube to the corresponding position in the output matrix
                            smallerY(i, j, k, :) = averaged_cube;
                        end
                    end
                end









                dvec = reshape(Y,[(size(Y,1)*size(Y,2)*size(Y,3)),size(Y,4)]);

                %remove rows of zeros
                dvec(~any(dvec,2),:) = [];
                dvec = double(dvec);

                nvec = size(dvec,1);

                dvec_dsampled(:,:) = interp1(1:nvec, dvec, linspace(1, nvec, nvoxels));

                % compute voxel-level connectivity matrix
                A(:,:,idx) = double(corr(dvec_dsampled(:,:)',dvec_dsampled(:,:)'));

                clear dvec data Y nvec dvec_dsampled dsY 
            end
        end
    end
end

save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/structured_input_data.mat');

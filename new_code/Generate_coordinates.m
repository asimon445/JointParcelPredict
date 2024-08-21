% This script is the first script that needs to be run to prepare the voxel
% level timeseries data for being fed into the functionally relevant
% parcellation algorithm
%
% It will identify the x, y, and z coordinates that cover by brain tissue 
% every subject. It will then remove voxels belonging to nodes in the
% cerebellum and subcortical structures. 

% Before running this script, be sure that there is a file that indicates
% where there is brain coverage in all subjects in the dataset (e.g., 
% 'HCP_overlap_mask.mat'). Be sure to also load in an atlas and a file
% indicating which nodes (if any) you want to remove from that atlas (e.g.,
% if you want to remove voxels from cerebellar or subcortical structures). 

clc; clear; close all;

addpath(genpath('/Users/ajsimon/Documents/Code/JointParcelPredict/'));

save_path = '/Users/ajsimon/Documents/Code/JointParcelPredict/';
remove_nodes = 0;   % set to 1 if we're removing nodes

Shen368 = load_untouch_nii('/Users/ajsimon/Documents/Code/JointParcelPredict/Shen_368_2mm.nii.gz');
Shen368 = Shen368.img;

coords_R = logical(false(91,109,91));
coords_L = logical(false(91,109,91));

load('/Users/ajsimon/Documents/Code/JointParcelPredict/HCP_overlap_mask.mat');
load('/Users/ajsimon/Documents/Code/JointParcelPredict/Shen368_10network.mat');

if remove_nodes == 1
    load('/Users/ajsimon/Documents/Code/JointParcelPredict/rmNodes_cortex.mat');
else
    rmNodes = [];
end

Shen368_nodeIx = [1:368]';

nodeIx = find(~ismember(Shen368_nodeIx(:,1),rmNodes)); 

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

clearvars -except coords_* save_path

% save the coordinates
save([save_path '/coordinates_to_subdivide_with_wholebrain.mat'],'coords_*');


clc; clear; close all;

files = dir('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Data/');
addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/sublist_and_bsiDep.mat');

% load in fMRI data and compute voxel-wise connectivity matrices
idx = 0;
for f = 1:length(files)
    isnii = strfind(files(f).name,'nii.gz');
    if ~isempty(isnii)
    
        for ff = 1:length(newsublist)

            submatch = strfind(lower(files(f).name),lower(newsublist{ff,1}));

            if ~isempty(submatch)
               
                idx=idx+1;

                if idx < 11
                    final_sublist{idx,1} = newsublist{ff,1};

                    d = load_untouch_nii([files(f).folder '/' files(f).name]);
                    data.img = d.img;

                    data.vec = reshape(data.img,[(size(data.img,1)*size(data.img,2)*size(data.img,3)),size(data.img,4)]);

                    %remove rows of zeros
                    data.vec(~any(data.vec,2),:) = [];
                    data.vec = double(data.vec);

                    outfile = sprintf('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/raw_timeseries/%s.mat',...
                        final_sublist{idx,1});

                    save(outfile,'data','-v7.3');
                    clear data
                end
            end
        end
    end
end
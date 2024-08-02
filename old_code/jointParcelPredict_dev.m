% This script will perform the joint prediction-parcellation 

clc; clear; close all;

files = dir('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Data/');
addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/sublist_and_bsiDep.mat');

nvoxels = 1000;   % This is the number of voxels we'll downsample to

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

save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/structured_input_data.mat')
%% 

% rho = 0.0001328;
rho = 0.00001;
K = 2;  % this is for creating the 2x2 piecewise constant

B = zeros(1000,1000);

check(1) = sum(B,'all');

for bi = 1:50

    [DL,L(bi,:)] = fGradient(squeeze(B(:,:,bi)),A,final_BSI_dep);  
    [res, B1, Z] = SpectralPart(DL, K);
    B(:,:,bi+1) = B(:,:,bi) - rho*res;

    clear DL

    check(bi+1) = -sum(B(:,:,bi+1),'all');
end




% try this with more people and real behavior. I should be able to see a
% convergence between zi and y at some point!

% for t=1:size(A,3)
%     %B(1) = 1;
%     B = ones(1000,1000);
%     for bi = 1:10
% %         zi(:,:,t) = A(:,:,t).*B(bi);   % zi should be 1 number for each patient (sum within individuals)
%         zi(:,:,t) = A(:,:,t).*B(:,:,bi);
% 
%         % compute loss function
% %         LBt = -(sum((beh(t)-zi(:,:,t))^2,"all"));
%         LBt = (beh(t)-sum(zi(:,:,t),'all'))^2;   
% 
%         dLBt = -2*zi(:,:,t)*A(:,:,t);
% 
%     end
% end
% 
% 






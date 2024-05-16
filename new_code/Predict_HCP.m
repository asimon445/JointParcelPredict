% This script will use the gradient boosting algorithm to predict IQ from 
% voxelwise connectivity matrices from the HCP dataset. 

clear; 

%% Setup the data for predictions
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

%% set parameters
K2 = 2; % K=2 for the first call
K3 = 3; % K=3 for the second call
rho = 0.5;
beta = 0.5;
iter_max = 100;
fixedStepSize = false;

nfolds = 2;
nperms = 1;

for p = 1:nperms

    indices = cvpartition(size(As,3),'k',nfolds);
                
    for f = 1:nfolds
        
        test.indx = indices.test(f);
        train.indx = indices.training(f);
        
        test.x = As(:,:,test.indx);
        train.x = As(:,:,train.indx);
        
        test.y = Y(indices.test(f),1);
        train.y = Y(indices.training(f),1);

        % Call SpectralGBR with K=2
        [testFit2_Fmat, testFit2_Fmat_list, testFit2_obj_seq, testFit2_EY_list, testFit2_rho_seq] = ...
            SpectralGBR(train.y, train.x, K2, rho, beta, iter_max, fixedStepSize);

        % Call SpectralGBR with K=3
        [testFit3_Fmat, testFit3_Fmat_list, testFit3_obj_seq, testFit3_EY_list, testFit3_rho_seq] = ...
            SpectralGBR(train.y, train.x, K3, rho, beta, iter_max, fixedStepSize);

        for fm = 1:length(testFit2_Fmat_list)
            EY_2(test.indx,fm,p) = AsByF(test.x,testFit2_Fmat_list{1,fm});
            EY_3(test.indx,fm,p) = AsByF(test.x,testFit3_Fmat_list{1,fm});
        end
    end
    
    % correlate estimated Y with actual Y
    for fm = 1:length(testFit2_Fmat_list)
        preds_k2(p,fm) = corr(Y,squeeze(EY_2(:,fm,p)));
        preds_k3(p,fm) = corr(Y,squeeze(EY_3(:,fm,p)));
    end

    if p == 1
        save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/pred_1perm_k2_100iter.mat','preds_k2');
        save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/pred_1perm_k3_100iter.mat','preds_k3');
    end
end

% save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/pred_10perms_k2_100iter.mat','preds_k2');
% save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/pred_10perms_k3_100iter.mat','preds_k3');
% 

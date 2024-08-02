% This will run regular CPM on the HCP data to establish a benchmark for
% prediction

clc; clear; close all;

%% set parameters and paths
addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

data_path = '/mnt/dustin/data/REST_LR/preprocessd_output/';
mask_path = '/mnt/dustin/data/REST_LR/matrix_new_vois/';

save_path = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/';

all_d = dir([data_path, '*_REST_LR_GSR_preprocessed_output.nii.gz']);
all_m = dir([mask_path, '*_LR_GSR_matrix_new_voi.nii.gz']);

no_sub = length(all_d);

addpath(genpath('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/Deviants/Deviants-selected/CPM/'));

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/HCP_IQ.mat');

nperms = 100;
nfolds = 10;

%% compute 268 node connectivity matrices
for d = 1:no_sub
    for m = 1:length(all_m)
        samesub = strcmp(all_d(d).name(1:6),all_m(m).name(1:6));
        if samesub == 1
            data = load_untouch_nii([all_d(d).folder, '/', all_d(d).name]);   % load data
            data = data.img;
            mask = load_untouch_nii([all_m(m).folder, '/', all_m(m).name]);   % load mask
            mask = mask.img;
            break;
        end
    end

    % turn data into 2D matrix -- dim 1 = location, dim 2 = time series
    dvec = reshape(data,[(size(data,1)*size(data,2)*size(data,3)),size(data,4)]);

    % turn mask into a vector
    mvec = reshape(mask,[],1);

    % compute connectivity matrices
    for n = 1:268   % loop through all nodes
        loc = find(mvec == n);
        ndata_TS(n,:) = mean(dvec(loc,:),1);
        clear loc
    end

    fprintf('%d ',d);
    tmpconn = corr(ndata_TS(:,:)',ndata_TS(:,:)');
    conns(:,:,d) = normalize(tmpconn);

    % find the Y data that matches this sub and format it for CPM
    for b = 1:length(HCP_IQ)
        IQsubnum = num2str(HCP_IQ(b,1));
        samesub = strcmp(IQsubnum,all_d(d).name(1:6));

        if samesub == 1
            Y(d,1) = HCP_IQ(b,2);
            break;
        end
        clear IQsubnum
    end
    clear ndata_TS tmpconn mvec dvec data mask

    if d == no_sub
        fprintf('\n');
    end
end

if length(Y) == size(conns,3)
    clearvars -except Y conns save_path nperms nfolds
    save([save_path '/HCP_268node_connmats_and_IQ.mat']);

    fprintf('Running predictions \n');

    for p = 1:nperms
        results = main(conns,Y,nfolds,'HCP');
        CPM_r(p,1) = results.r_rank;
        clear results

        fprintf('%d ',p);
        if p == nperms
            fprintf('\n');
        end
    end

    save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Test_results/CPM_results.mat','CPM_r');

else
    error('The number of connectivity matrices does not match the number of subjects with IQ data');
end


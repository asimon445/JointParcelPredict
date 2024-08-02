% This script will feed the MxNxP voxelwise connectivity matrices and the
% vector of behavioral data into the gradient descent and spectral
% clustering algorithms. With each iteration, we will seek to find the
% optimal step size for minimizing the loss. 

% Before running this script, you will first need to run 'setup_data.m'

% parameters
rho = 1;
K = 2;  % this is for creating the 2x2 piecewise constant
niter = 50;  % specify the number of iterations 

% load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/structured_input_data.mat');

ix=0;
fs = 0;  % set this to zero when an optimal step size has been identified to exit the while loop

while fs == 0
    ix=ix+1;
    B(:,:) = zeros(1000,1000);

    loss(ix,1) = sum(B(:,:),'all');

    for bi = 1:niter

        [DL,L(bi,:)] = fGradient(squeeze(B(:,:,bi)),A,final_BSI_dep);
        [res, B1, Z] = SpectralPart(DL, K);
        B(:,:,bi+1) = B(:,:,bi) - rho*res;

        clear DL

        loss(ix,bi+1) = -sum(B(:,:,bi+1),'all');
    end

    % determine if loss is minimized at some point during the iterations
    fminsearch(y,1);



%     rho_new = backtr(rho,ix,1,loss);

end


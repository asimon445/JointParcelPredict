% this is just code for testing out the joint parcellation-prediction stuff

% to do:
% 1. double check that downsampling is taking the average of the
% neighboring voxels, and NOT just deleting stuff to make the matrix its
% desired size
%
% 2. Do we need 'Z' from 'SpectralPart'?



clear; 

path = '/Users/ajsimon/Documents/Data/Constable_lab/Flexible_parcellations_CPM/dev_data/raw/';

%load('/Users/ajsimon/Documents/Data/Constable_lab/Flexible_parcellations_CPM/dev_data/B.mat');
load('/Users/ajsimon/Documents/Data/Constable_lab/Flexible_parcellations_CPM/dev_data/beh.mat');

files = dir(path);

nvoxels = 1000;   % This is the number of voxels we'll downsample to
ds_fac = 3;   % downsample each axis by a factor of this number  - will reduce the matrix size by n^3

idx = 0;
for f = 1:length(files)
    if files(f).isdir == 0
        idx=idx+1;
        data = load_untouch_nii([path files(f).name]);

        Y = data.img;

%         % get rid of first row/col of Y for x,y,z coords, since it's just
%         % 0's anyway and it makes downsampling tricky
%         Y(1,:,:,:) = [];
%         Y(:,1,:,:) = [];
%         Y(:,:,1,:) = [];
% 
%         % downsample errthang 
%         d = size(Y);
%         
%         for s = 1:size(Y,4)
%             X = reshape(Y(:,:,:,s), ds_fac, d(1)/ds_fac, ds_fac, d(2)/ds_fac, ds_fac, d(3)/ds_fac);
%             
%             dsY(:,:,:,s) = squeeze(sum(X, [1,3,5])) / 8;
%             clear X
%         end
% 
%         % transform voxel level time series into a vector (i guess matrix, where
%         % dim 2 is time)
%         dvec = reshape(dsY,[(size(dsY,1)*size(dsY,2)*size(dsY,3)),size(dsY,4)]);

        dvec = reshape(Y,[(size(Y,1)*size(Y,2)*size(Y,3)),size(Y,4)]);
        
        %remove rows of zeros
        dvec(~any(dvec,2),:) = [];  
        dvec = double(dvec);

        nvec = size(dvec,1);

        dvec_dsampled(:,:) = interp1(1:nvec, dvec, linspace(1, nvec, nvoxels));

        % compute voxel-level connectivity matrix
%         eval(sprintf('A.num%d(:,:) = double(corr(dvec(:,:)'',dvec(:,:)''));',idx));   %dvec needs to be reduced by ~10x for matlab to be able to doo this  

        A(:,:,idx) = double(corr(dvec_dsampled(:,:)',dvec_dsampled(:,:)'));
        % compute dummy behavioral data
%         beh(idx,1) = rand(1,1);

        clear dvec data Y nvec dvec_dsampled dsY
    end
end


% % Compute a matrix of the same dimensions as A that contains the prediction
% % value of each voxel -- it will eventually be a prediction, but for now
% % just correlate each edge with behavior just for development purposes
% for r = 1:nvoxels
%     B(r,:) = corr(squeeze(A(r,:,:))',beh);
% end
% 
% 
% % replace all values on the diagonal with a 0 -- this is a temp solution to
% % replacing the NaNs on the diagonal
% for d = 1:15000
%     B(d,d) = 0;
% end

% start with B = 0

% compute mean of A
Am = mean(A,3);

loss0 = 0;

% z = A*B for each sub
rho = 1;
K = 2;  % this is for creating the 2x2 piecewise constant

for s=1:size(A,3)
    ix=0;
    for b = 0:100:10000
        ix=ix+1;
        temp = A(:,:,s).*b';   % zi should be 1 number for each patient (sum within individuals)
        zi = sum(temp,"all");
        clear temp

        % compute loss function
        loss(s) = loss0-(beh(s) - zi)^2;
        dloss(:,:,ix,s) = -2*zi*A(:,:,s);

        % do the spectal part
        [res(:,:,ix,s), B1(:,:,ix,s), Z] = SpectralPart(squeeze(dloss(:,:,ix,s)), K);
                
        Bnew(:,:,ix,s) = b+B1(:,:,ix,s);

        test1(ix) = Bnew(1,1,ix,s);
        test2(ix) = Bnew(1,2,ix,s);
        test3(ix) = Bnew(2,1,ix,s);
        test4(ix) = Bnew(2,2,ix,s);

        clear zi Z
    end
end














% clear
% clc
% 
% load('sample_sparse');
% 
% 
% index_nii=load_nii_gz('sample_indexmap.nii.gz', '/mridata2/mri_group/xilin_data/xenios/');
% ind_img = index_nii.img;
% clear index_nii;
% 
% nozero = find(ind_img);
% len = length(nozero);
% 
% w = sparse(sparseM(:,1), sparseM(:,2), exp(- (sparseM(:,3)/2).^2), len, len);
% 
% no_cluster = 50;
% 
% [Eigenvectors,Eigenvalues, vbar] = my_ncut(w,no_cluster);
% 
% save('anatomical_eig', 'Eigenvectors');

% compute discretize Ncut vectors
[Discrete, oEigenvectors, ini_label, ncvalue] =discretisation(Eigenvectors(:, 1:20));


save('anatomical_discre', 'Discrete');

dim = size(ind_img);
label_img = zeros(dim);

for i=1:len;
    cur = full(Discrete(i,:));
    aa = find(cur);
    if( length(aa)==1)
        label_img(nozero(i)) = aa;
    else
        display(['something wrong for voxel ', num2str(i)]);
    end
end

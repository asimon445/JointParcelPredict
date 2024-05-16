

clc; clear; close all;

addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

data_path = '/mnt/dustin/data/REST_LR/preprocessd_output/';
save_path = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/';

all_d = dir([data_path, '*_REST_LR_GSR_preprocessed_output.nii.gz']);

no_sub = length(all_d);

for i = 1: no_sub
    nii = load_untouch_nii([all_d(i).folder, '/', all_d(i).name]); 
    img = nii.img;
    img_norm = sum(img.^2, 4);
    
    if( i==1)
        dim = size(img_norm);
        mask_all = zeros([dim, no_sub]);
    end
    mask_all(:,:,:,i)= double(img_norm>0);
end

mask_all= uint8(mask_all);

save([save_path, 'HCP_overlap_mask'], 'mask_all');


% dpath = '/mnt/dustin/data/REST_LR/preprocessd_output/';
% % dpath = '/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Data/';
% 
% files = dir(dpath);
% 
% %% mask out voxels that are not present in every single subject
% ix=0;
% for f = 1:length(files)
% 
%     isvox = strfind(files(f).name,'REST_LR_GSR_preprocessed_output.nii.gz');
% %     isvox = strfind(files(f).name,'bold.nii.gz');
%     
%     if ~isempty(isvox)
%         ix=ix+1;
%         file = sprintf('%s%s',dpath,files(f).name);
%         data = load_untouch_nii(file);
%         sublist{ix,1} = files(f).name(1:6);
% 
%         Y = data.img;
% 
%         % take the sum, because if ALL values in the timeseries are zero then
%         % that voxel must not contain brain data
%         image_sums(:,:,:,ix) = sum(Y,4);
% 
%         clear vec Y data file
% 
%         fprintf('file number %d \n',ix);
%         if ix == 100
%             break;
%         end
%     end
% end
% 
% save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/HCP_LR_indiv_voxel_coverage_mask.mat','image_sums','-v7.3');
% 
non_zero_coords = true(size(image_sums,1:3));

for i = 1:size(image_sums, 1)
    for j = 1:size(image_sums, 2)
        for k = 1:size(image_sums, 3)
            % Check if the pixel is non-zero across all images
            if any(image_sums(i,j,k,:) == 0)
                non_zero_coords(i,j,k) = false;
            end
        end
    end
end
% 
% % Index the non-zero coordinates
% [x_pos, y_pos, z_pos] = ind2sub(size(image_sums,1:3), find(non_zero_coords));

%% 
ix=0;
for f = 1:length(all_d)

    isvox = strfind(all_d(f).name,'REST_LR_GSR_preprocessed_output.nii.gz');

    if ~isempty(isvox)
        ix=ix+1;
        file = sprintf('%s%s',data_path,all_d(f).name);
        data = load_untouch_nii(file);

        img = data.img;

        % populate this patient's image with only voxels that everyone has
        img_cleaned = zeros(size(img));
        for v = 1:size(x_pos)
            img_cleaned(x_pos(v),y_pos(v),z_pos(v),:) = img(x_pos(v),y_pos(v),z_pos(v),:);
        end

        % average together voxels in a 4x4x4 lattice to reduce the size of the full
        % matrix
        ixx=0;
        for x = 2:4:size(img_cleaned,1)-2
            ixx=ixx+1;
            ixy=0;
            for y = 2:4:size(img_cleaned,2)-2
                ixy=ixy+1;
                ixz=0;
                for z = 2:4:size(img_cleaned,3)-2
                    ixz=ixz+1;

                    centerTS = squeeze(img_cleaned(x,y,z,:))';
                    cube = img_cleaned(x-1:x+2,y-1:y+2,z-1:z+2,:);

                    % take the neighboring average if the center of the cube has
                    % timeseries data that isn't all zeros
                    numvals = length(unique(centerTS));
                    if numvals > 1

                        % replace empty voxels with nans
                        for xx=1:4
                            for yy=1:4
                                for zz=1:4
                                    temp = squeeze(cube(xx,yy,zz,:));
                                    temp(~any(temp,2),:) = NaN;
                                    cube(xx,yy,zz,:) = temp;
                                    clear temp
                                end
                            end
                        end
                        Y_reduced(ixx,ixy,ixz,:) = mean(cube,[1,2,3],'omitnan');

                    else
                        Y_reduced(ixx,ixy,ixz,1:size(img_cleaned,4)) = 0;
                    end

                    clear centerTS cube numvals
                end
            end
        end

        dvec = reshape(Y_reduced,[(size(Y_reduced,1)*size(Y_reduced,2)*size(Y_reduced,3)),size(Y_reduced,4)]);

        %remove rows of zeros
        dvec(~any(dvec,2),:) = [];
        dvec = double(dvec);

        connmat(:,:) = double(corr(dvec(:,:)',dvec(:,:)'));

        outfile = sprintf('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/voxelwise_connmats/HCP_%s_4vox.mat',...
            all_d(ix).name(1:11));

        save(outfile,'connmat');
        clear dvec Y_reduced connmat outfile

        % average together voxels in a 3x3x3 lattice to reduce the size of the full
        % matrix
        ixx=0;
        for x = 2:3:size(img_cleaned,1)-1
            ixx=ixx+1;
            ixy=0;
            for y = 2:3:size(img_cleaned,2)-1
                ixy=ixy+1;
                ixz=0;
                for z = 2:3:size(img_cleaned,3)-1
                    ixz=ixz+1;

                    centerTS = squeeze(img_cleaned(x,y,z,:))';
                    cube = img_cleaned(x-1:x+2,y-1:y+2,z-1:z+2,:);

                    % take the neighboring average if the center of the cube has
                    % timeseries data that isn't all zeros
                    numvals = length(unique(centerTS));
                    if numvals > 1

                        % replace empty voxels with nans
                        for xx=1:3
                            for yy=1:3
                                for zz=1:3
                                    temp = squeeze(cube(xx,yy,zz,:));
                                    temp(~any(temp,2),:) = NaN;
                                    cube(xx,yy,zz,:) = temp;
                                    clear temp
                                end
                            end
                        end
                        Y_reduced(ixx,ixy,ixz,:) = mean(cube,[1,2,3],'omitnan');

                    else
                        Y_reduced(ixx,ixy,ixz,1:size(img_cleaned,4)) = 0;
                    end

                    clear centerTS cube numvals
                end
            end
        end

        dvec = reshape(Y_reduced,[(size(Y_reduced,1)*size(Y_reduced,2)*size(Y_reduced,3)),size(Y_reduced,4)]);

        %remove rows of zeros
        dvec(~any(dvec,2),:) = [];
        dvec = double(dvec);

        connmat(:,:) = double(corr(dvec(:,:)',dvec(:,:)'));

        outfile = sprintf('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/voxelwise_connmats/HCP_%s_3vox.mat',...
            all_d(ix).name(1:11));

        save(outfile,'connmat');
        clear dvec Y_reduced connmat outfile
    end
end

% save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/connmats.mat','connmat','-v7.3');
% save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/connmats_sublist.mat','sublist');

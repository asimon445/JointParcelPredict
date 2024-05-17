% This script will generate matrices containing information about where
% the voxels are in space. 

nii = load_untouch_nii('/Users/ajsimon/Documents/Data/Constable_lab/Flexible_parcellations_CPM/dev_data/old/raw/test_data1.nii.gz');
img = nii.img;

load('/Users/ajsimon/Documents/Data/Constable_lab/Flexible_parcellations_CPM/dev_data/data/Get_spatial/HCP_overlap_mask.mat');

img_cleaned = zeros(size(img));
for v = 1:size(x_pos)
    img_cleaned(x_pos(v),y_pos(v),z_pos(v),:) = img(x_pos(v),y_pos(v),z_pos(v),:);
end

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


img_norm = sum(Y_reduced.^2, 4);

nzero_reduced = true(size(img_norm,1:3));

img_norm = sum(Y_reduced.^2, 4);
mask_reduced = double(img_norm >0);
nzero_reduced = true(size(img_norm,1:3));
for i = 1:size(img_norm, 1)
    for j = 1:size(img_norm, 2)
        for k = 1:size(img_norm, 3)
            % Check if the pixel is non-zero across all images
            if any(img_norm(i,j,k) == 0)
                nzero_reduced(i,j,k) = false;
            end
        end
    end
end

% Index the non-zero coordinates
[x_red, y_red, z_red] = ind2sub(size(img_norm,1:3), find(nzero_reduced));

% Reshape Y_reduced to a 2D matrix
dvec = reshape(Y_reduced, [(size(Y_reduced, 1) * size(Y_reduced, 2) * size(Y_reduced, 3)), size(Y_reduced, 4)]);

% Remove rows that are all zeros
dvec(~any(dvec, 2), :) = [];
dvec = double(dvec);

% 
num_non_zero = length(x_red);
neighbormat = zeros(num_non_zero, num_non_zero);

% Step 8: Check for neighbors
for i = 1:num_non_zero
    for j = 1:num_non_zero
        if i ~= j
            % Check if positions are neighbors (26 possibilities)
            if abs(x_red(i) - x_red(j)) == 1 && abs(y_red(i) - y_red(j)) == 0 && abs(z_red(i) - z_red(j)) == 0  
                neighbormat(i, j) = 1;
            elseif abs(x_red(i) - x_red(j)) == 0 && abs(y_red(i) - y_red(j)) == 1 && abs(z_red(i) - z_red(j)) == 0  
                neighbormat(i, j) = 1;
            elseif abs(x_red(i) - x_red(j)) == 0 && abs(y_red(i) - y_red(j)) == 0 && abs(z_red(i) - z_red(j)) == 1  
                neighbormat(i, j) = 1;
            elseif abs(x_red(i) - x_red(j)) == 1 && abs(y_red(i) - y_red(j)) == 1 && abs(z_red(i) - z_red(j)) == 0
                neighbormat(i, j) = 1;
            elseif abs(x_red(i) - x_red(j)) == 1 && abs(y_red(i) - y_red(j)) == 0 && abs(z_red(i) - z_red(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_red(i) - x_red(j)) == 0 && abs(y_red(i) - y_red(j)) == 1 && abs(z_red(i) - z_red(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_red(i) - x_red(j)) == 1 && abs(y_red(i) - y_red(j)) == 1 && abs(z_red(i) - z_red(j)) == 1
                neighbormat(i, j) = 1;
            end
        end
    end
end

% visualize neighbormat
figure; imagesc(neighbormat);

% create a vector containing info about total number of neighbors each
% voxel has
total_neighbors = sum(neighbormat, 2);
figure; scatter([1:1644],total_neighbors)




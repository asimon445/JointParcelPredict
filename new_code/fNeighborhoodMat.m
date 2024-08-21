function [neighbormat,total_neighbors] = fNeighborhoodMat(subvolume_coords)
% This function will compute a neighborhood matrix and a vector containing
% info about how many neighboring voxels each voxel contains to be used in
% the functional parcellation algorithm.

x_pos = subvolume_coords(:,1);
y_pos = subvolume_coords(:,2);
z_pos = subvolume_coords(:,3);

num_non_zero = length(x_pos);
neighbormat = zeros(num_non_zero,num_non_zero);

for i = 1:num_non_zero
    for j = 1:num_non_zero
        if i ~= j
            % Check if positions are neighbors (26 possibilities)
            if abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 0 && abs(z_pos(i) - z_pos(j)) == 0
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 0 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 0
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 0 && abs(y_pos(i) - y_pos(j)) == 0 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 0
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 0 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 0 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            elseif abs(x_pos(i) - x_pos(j)) == 1 && abs(y_pos(i) - y_pos(j)) == 1 && abs(z_pos(i) - z_pos(j)) == 1
                neighbormat(i, j) = 1;
            end
        end
    end
end

total_neighbors = sum(neighbormat, 2);

end
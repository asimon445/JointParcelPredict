function [expanded_subvolumes] = fExpandSubvolumes(subvolumes, mask)
    % Initialize the cell array to store the expanded subvolumes
    expanded_subvolumes = cell(size(subvolumes));

    % Get the size of the mask to ensure we don't exceed boundaries
    [mask_x, mask_y, mask_z] = size(mask);

    % Define the neighborhood for a 1-voxel thick expansion
    [X, Y, Z] = ndgrid(-1:1, -1:1, -1:1);
    neighbor_offsets = [X(:) Y(:) Z(:)];

    for i = 1:length(subvolumes)
        current_subvolume = subvolumes{i};
        expanded_voxels = [];

        % For each voxel in the subvolume, add the neighboring voxels
        for j = 1:size(current_subvolume, 1)
            voxel = current_subvolume(j, :);
            neighbors = bsxfun(@plus, voxel, neighbor_offsets);

            % Ensure neighbors are within mask boundaries
            valid_neighbors = neighbors(all(neighbors > 0, 2) & ...
                                        neighbors(:, 1) <= mask_x & ...
                                        neighbors(:, 2) <= mask_y & ...
                                        neighbors(:, 3) <= mask_z, :);

            % Check if the neighbors are within the brain mask
            mask_indices = sub2ind(size(mask), valid_neighbors(:,1), ...
                                   valid_neighbors(:,2), valid_neighbors(:,3));
            valid_neighbors = valid_neighbors(mask(mask_indices) == 1, :);

            % Add the valid neighbors to the expanded_voxels list
            expanded_voxels = [expanded_voxels; valid_neighbors]; %#ok<AGROW>
        end

        % Remove duplicates by converting to unique rows
        expanded_voxels = unique(expanded_voxels, 'rows');
        
        % Store the expanded subvolume
        expanded_subvolumes{i} = expanded_voxels;
    end
end
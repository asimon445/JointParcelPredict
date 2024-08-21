function [new_subvolumes] = fMergeSubvolumes(subvolumes, min_voxels,max_voxels)
    % This script will merge subvolumes of the brain based off euclidean 
    % distances. Merging will be done only on subvolumes that contain fewer
    % than a specified number of voxels ('min_voxels')
    % 
    % Input arguments:
    %  'subvolumes': a cell array of x, y, and z coordinates for all
    %       subvolumes of the brain
    %  'min_voxels': all subvolumes with fewer voxels than this number will
    %       get merged with a neighboring subvolume
    %  'max_voxels': any subvolume with more voxels than this number will
    %       have small subvolumes merged with it.

    nBordering = zeros(length(subvolumes)); % Initialize nBordering matrix

    for i = 1:length(subvolumes)
        current_idx = i;
        current_subvolume = subvolumes{i};

        % Loop through all other subvolumes to find the number of bordering
        % voxels between all possible pairs of subvolumes
        for j = 1:length(subvolumes)
            if j ~= i  
                other_subvolume = subvolumes{j};
                distances = pdist2(current_subvolume, other_subvolume, 'euclidean');
                nBordering(i,j) = length(find(distances == 1));
            end
        end

        % find the right subvolume to merge to
        %max_nBordering(i) = max(nBordering(i,:));
        [a,b] = sort(nBordering(i,:),'descend');
        b(a==0) = [];
        a(a==0) = [];

        found = 0; inix = 0;
        while found == 0
            try
                inix=inix+1;
                if length(subvolumes{b(inix)}) < max_voxels
                    posMaxBordering(i) = b(inix);
                    found = 1;
                end
            catch
                error('All neighboring subvolumes contain more than the max number of voxels to merge into');
            end
        end
        clear a b
    end

    % sort by size and handle the really small subvolumes first
    lens = cellfun(@length, subvolumes);

    % find all subvolumes that are small and need to be merged
    smallSubvols = find(lens < min_voxels); 

    % list the nearest subvolumes to the ones that need to be merged
    mergeTo = posMaxBordering(smallSubvols);   

    % find the ones that exist in mergeTo and smallSubvols. Remove it from
    % smallSubvols.
    ix_inboth = intersect(smallSubvols, mergeTo);
    merged_inBoth = [];
    new_subvolumes = {};

    for bth = 1:length(ix_inboth)
        inBoth_smallSubvol_ix = find(smallSubvols == ix_inboth(bth));
        inBoth_mergeTo_ix = find(mergeTo == ix_inboth(bth));

        smallSubvol_mergeSmall_ix = smallSubvols(inBoth_mergeTo_ix);
        smallSubvol_rep_ix = smallSubvols(inBoth_smallSubvol_ix);

        if ~isempty(smallSubvol_rep_ix) && ~isempty(smallSubvol_mergeSmall_ix)
            new_subvolumes{end+1} = [subvolumes{smallSubvol_mergeSmall_ix}; subvolumes{smallSubvol_rep_ix}];
            merged_inBoth(end+1:end+2) = [smallSubvol_mergeSmall_ix, smallSubvol_rep_ix];

            % Remove the merged ones from both mergeTo and smallSubvol
            mergeTo([inBoth_mergeTo_ix, inBoth_smallSubvol_ix]) = [];
            smallSubvols([inBoth_mergeTo_ix, inBoth_smallSubvol_ix]) = [];
        end
    end

    % merge unique ones first
    [uniqueMergeTos, ~, uniqueMergeTo_ixs] = unique(mergeTo, 'stable');

    for s = 1:length(uniqueMergeTos)
        ixs = find(mergeTo == uniqueMergeTos(s));
        combined_subvolume = subvolumes{uniqueMergeTos(s)};
        
        for i = 1:length(ixs)
            smallSubvols_ix = smallSubvols(ixs(i));
            combined_subvolume = [combined_subvolume; subvolumes{smallSubvols_ix}];
        end
        
        new_subvolumes{end+1} = combined_subvolume;
    end

    % identify the subvolumes that did not have anything merged into them
    for s = 1:length(subvolumes)
        issmall = find(smallSubvols == s);
        isMergedTo = find(mergeTo == s);
        isMerged_inBoth = find(merged_inBoth == s);

        if isempty(issmall) && isempty(isMergedTo) && isempty(isMerged_inBoth)
            new_subvolumes{end+1} = subvolumes{s};
        end
    end

end
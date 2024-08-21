% This script will prepare the voxel time series data to be inputted into
% the joint parcellation-prediction algorithm. 

clear; 

load('/Users/ajsimon/Documents/Code/JointParcelPredict/coordinates_to_subdivide_with_incl_subcort_exl_cereb.mat');

side = 'L';
outdir = sprintf('/Users/ajsimon/Documents/Data/Constable_lab/Functional_parcellations/HCP/Subvolumes_%s/',side);
outdir_n = sprintf('/Users/ajsimon/Documents/Data/Constable_lab/Functional_parcellations/HCP/Subvolume_neighborhoods_%s/',side);

coords = eval(sprintf('coords_%s',side));
[x_pos, y_pos, z_pos] = ind2sub(size(coords,1:3), find(coords));

% find min/max non-zero coords in each direction
min_x = min(x_pos);   % min/max x coords will differ between hemispheres
max_x = max(x_pos);
min_y = min(y_pos);
max_y = max(y_pos);
min_z = min(z_pos);
max_z = max(z_pos);

inc_x = round((max_x - min_x) / 4);
inc_y = round((max_y - min_y) / 4);
inc_z = round((max_z - min_z) / 3);

% Initialize a cell array to store the subvolumes
subvolumes = {}; % Loop through the matrix with the defined step size

% break the matrix into smol cubes
for z = min_z:inc_z:max_z   % superior/inferior
    for y = min_y:inc_y:max_y   % front/back
        for x = min_x:inc_x:max_x   % Left/right

            zrng = find(z_pos >= z & z_pos < z+inc_z);
            yrng = find(y_pos >= y & y_pos < y+inc_y);
            xrng = find(x_pos >= x & x_pos < x+inc_x);

            % Find the common numbers
            sharedCoords = intersect(intersect(xrng, yrng), zrng);

            if ~isempty(sharedCoords)

                for i = 1:length(sharedCoords)
                    subvolume(i,1) = x_pos(sharedCoords(i,1));
                    subvolume(i,2) = y_pos(sharedCoords(i,1));
                    subvolume(i,3) = z_pos(sharedCoords(i,1));
                end
                
                subvolumes{end+1} = subvolume;
                clear subvolume
            end
        end
    end
end

%% Combine subvolumes iteratively -- you'll need to play with these parameters. These seem to work with HCP data
if contains(side,'L')
    new_subvolumes = fMergeSubvolumes(subvolumes,700,2500);
    new_subvolumes = fMergeSubvolumes(new_subvolumes,1000,1800);
    new_subvolumes = fMergeSubvolumes(new_subvolumes,1300,1800);
elseif contains(side,'R')
    new_subvolumes = fMergeSubvolumes(subvolumes,700,2500);
    new_subvolumes = fMergeSubvolumes(new_subvolumes,800,2000);
    new_subvolumes = fMergeSubvolumes(new_subvolumes,1200,1800);
end

% check that the subvolume merging was done correctly 
nvox = length(find(coords == 1));

for i = 1:length(new_subvolumes)
    newlens(i) = length(new_subvolumes{i});
end

if sum(newlens) ~= nvox
    error('Subvolume merging was done incorrectly');
end

% Expand each subvolume so that it overlaps a little with neighboring subvolumes
expanded_subvolumes = fExpandSubvolumes(new_subvolumes, coords);

% Compute neighborhood matrices for each subvolume and save
for s = 1:length(expanded_subvolumes)

    subvol = expanded_subvolumes{s};
    [neighbormat,total_neighbors] = fNeighborhoodMat(subvol);
    
    outfile = sprintf('%sSubvol_coords_%d.mat',outdir,s);
    outfile_n = sprintf('%sSubvol_neighborhoodInfo_%d.mat',outdir_n,s);

    save(outfile,'subvol');
    save(outfile_n,'neighbormat','total_neighbors');
    
    clear subvol neighbormat total_neighbors outfile*
end


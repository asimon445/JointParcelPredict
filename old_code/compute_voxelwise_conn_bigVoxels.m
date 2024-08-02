clc; clear; close all;

files = dir('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/Data/');
addpath(genpath('/home/ajs332/Documents/NIfTI_20140122/'));

load('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/sublist_and_bsiDep.mat');

% load in fMRI data and compute voxel-wise connectivity matrices
idx = 0;
for f = 1:length(files)
    isnii = strfind(files(f).name,'nii.gz');
    if ~isempty(isnii)
    
        for ff = 1:length(newsublist)

            submatch = strfind(lower(files(f).name),lower(newsublist{ff,1}));

            if ~isempty(submatch)
                idx=idx+1;
                final_sublist{idx,1} = newsublist{ff,1};
                final_BSI_dep(idx,1) = new_bsi(ff,1);

                data = load_untouch_nii([files(f).folder '/' files(f).name]);
                Y = data.img;

                % average together voxels in a 3x3x3 lattice to reduce the size of the full
                % matrix
                ixx=0;
                for x = 2:3:size(Y,1)-1
                    ixx=ixx+1;
                    ixy=0;
                    for y = 2:3:size(Y,2)-1
                        ixy=ixy+1;
                        ixz=0;
                        for z = 2:3:size(Y,3)-1
                            ixz=ixz+1;

                            centerTS = squeeze(Y(x,y,z,:))';

                            cube = Y(x-1:x+1,y-1:y+1,z-1:z+1,:);

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
                                Y_reduced(ixx,ixy,ixz,1:size(Y,4)) = 0;
                            end

                            clear centerTS cube numvals
                        end
                    end
                end

                dvec = reshape(Y_reduced,[(size(Y_reduced,1)*size(Y_reduced,2)*size(Y_reduced,3)),size(Y_reduced,4)]);

                %remove rows of zeros
                dvec(~any(dvec,2),:) = [];
                dvec = double(dvec);

                connmat = double(corr(dvec(:,:)'',dvec(:,:)''));

                outfile = sprintf('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/voxelwise_connmats/%s.mat',...
                    final_sublist{idx,1});

                save(outfile,'connmat');

                clear dvec Y Y_reduced data connmat outfile
                fprintf('%d ',idx);
            end
        end
    end
end

save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/BSI_dep.mat','final_BSI_dep');
save('/data22/mri_group/dustinlab_data/dustinlab/Documents/AJ/JointParcelPredict_dev/sublist.mat','final_sublist');



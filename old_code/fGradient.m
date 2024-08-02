function [DL,L] = fGradient(B,A,y)
% Input arguments:
%  B -- an nxm matrix of zeros (for first iteration), and results from the
%       spectral part + B for every other iteration
%  A -- an nxmxp matrix of voxelwise connectivity for all patients
%  y -- a vector of behavioral data for all patients
% 

for t = 1:size(A,3)

    zi = squeeze(A(:,:,t)).*B(:,:);

    z = sum(zi,'all');

    tempL(t) = (y(t) - z)^2;
    tempDL(:,:,t) = (z-y(t))*squeeze(A(:,:,t));
%     tempDL(:,:,t) = (y(t)-z)*squeeze(A(:,:,t));
    clear zi z

end

L = sum(tempL);
DL = 2*sum(tempDL,3);

end


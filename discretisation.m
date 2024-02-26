function [EigenvectorsDiscrete,EigenVectors, ini_label, lastObjectiveValue]=discretisation(EigenVectors)
% 
% EigenvectorsDiscrete=discretisation(EigenVectors)
% 
% Input: EigenVectors = continuous Ncut vector, size = ndata x nbEigenvectors 
% Output EigenvectorsDiscrete = discrete Ncut vector, size = ndata x nbEigenvectors
%
% Timothee Cour, Stella Yu, Jianbo Shi, 2004

[n,k]=size(EigenVectors);

vm = sqrt(sum(EigenVectors.*EigenVectors,2));
EigenVectors = EigenVectors./repmat(vm,1,k);


% generate the R matrix by random initialization
R=zeros(k);
R(:,1)=EigenVectors(1+round(rand(1)*(n-1)),:)';
c=zeros(n,1);
for j=2:k
    c=c+abs(EigenVectors*R(:,j-1));
    [minimum,i]=min(c);
    R(:,j)=EigenVectors(i,:)';
end
% the number of column in R is the same as the number of clusters

lastObjectiveValue=0;
exitLoop=0;
nbIterationsDiscretisation = 0;
nbIterationsDiscretisationMax = 20;%voir

ini = 0;

while exitLoop== 0 
    nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;
    % assign the labels by maximum suppression
    EigenvectorsDiscrete = discretisationEigenVectorData(EigenVectors*R);  % binary matrix
    
    % this part was added to the code to get the initial labeling 
    if(ini ==0)
        ini_label = zeros(1,n);
        for i=1:n;
            ini_label(i) = find(EigenvectorsDiscrete(i,:));
        end
        ini=1;
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % calculate the R matrix
    [U,S,V] = svd(EigenvectorsDiscrete'*EigenVectors,0);    
    NcutValue=2*(n-trace(S));
    
    if abs(NcutValue-lastObjectiveValue) < eps | nbIterationsDiscretisation > nbIterationsDiscretisationMax
        exitLoop=1;
    else
        lastObjectiveValue = NcutValue;
        R=V*U';
    end
end
% Function to calculate the As array inner product with F 
function res = AsByF(As, F)
% Extracting the size of the As array and F matrix
[p, ~, N] = size(As);

% Initializing the result vector
res = NaN(1, N);

% Looping through the third dimension of As
for i = 1:N
    % Multiplying element-wise and summing
    res(i) = sum(sum(As(:,:,i) .* F));
end
end
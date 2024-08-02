

function [res, B, Z] = SpectralPart(Y, K)
   % [row, col] = find(isnan(YourMatrix));    
    [U, S, V] = svds(Y, K);
    [cluster, ~] = kmeans(V, K, 'MaxIter', 50, 'Replicates', 30);
    n = size(Y, 1);
    Z = sparse(1:n, cluster, 1, n, K);
    B = zeros(K, K);
    for j = 1:K
        for k = j:K
            mask_j = cluster == j;
            mask_k = cluster == k;
            mask = mask_j & mask_k';

            if any(mask(:))
                B(j, k) = mean(Y(mask));
            else
                B(j, k) = 0; % Assign a default value if no elements are selected
            end
            B(k, j) = B(j, k);
        end
    end
    ZB = Z * B;
    res = ZB * Z';
end





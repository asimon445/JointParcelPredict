% Function to calculate the loss and the gradient
function [loss, grad] = LossAndGradient(Y, As, F)
    % Extracting the number of observations
    N = length(Y);
    
    % Initializing loss and gradient
    loss = 0;
    [p, ~] = size(F);
    grad = zeros(p, p);
    
    % Looping through each observation
    for i = 1:N
        % Calculating the difference for each observation
        tmp_diff = sum(sum(As(:,:,i) .* F)) - Y(i);
        
        % Updating loss
        loss = loss + tmp_diff^2;
        
        % Updating gradient
        grad = grad + 2 * tmp_diff * As(:,:,i);
    end
end

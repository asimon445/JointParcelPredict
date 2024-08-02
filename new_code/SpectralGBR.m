% As is the array of all adjacency matrices. The 3rd dimension is the index of observation
% Y is the vector of all response variables.
% K is the number of partitions in each update (the basic version is K=2. Usually, K>2 usually improves the performance, but it may be slower. It depends on the scale of the problem.)
% rho is the step size of the gradient step
% beta is the shrinkage factor for backtracking line search
% We are solving the problem 
% min_{F} sum_{i} \sum_{i=1}^N (Y[i] - sum(A[,,i]*F))^2
% s.t. F is a symmetric blockwise constant matrix
% The algorithm is a gradient boosting algorithm, rather than projected gradient descent (another option)




function [Fmat, Fmat_list, obj_seq, EY_list, rho_seq] = SpectralGBR(Y, As, K, rho, beta, iter_max, fixedStepSize, lambda, LR)
    if nargin < 7
        fixedStepSize = false; % Default value if not specified
    end
    if nargin < 6
        iter_max = 50; % Default value if not specified
    end
    
    N = length(Y);
    p = size(As, 1);
    Fmat = zeros(p, p);
    [loss, grad] = LossAndGradient(Y, As, Fmat, lambda, LR);
    obj_seq = inf(1, iter_max);
    obj_seq(1) = loss;
    rho_seq = zeros(1, iter_max);
    Fmat_list = cell(1, iter_max);
    EY_list = cell(1, iter_max);
    
    for iter = 1:iter_max
        fprintf('iteration %d\n', iter);
        fprintf('%f\n', loss);
        Fmat_old = Fmat;
        step_size_OK = false;
        
        % Keep reducing the step size until the loss is reduced or the step size is too small
        while ~step_size_OK && rho > 1e-8
            projected_grad = SpectralPart(grad, K); % Needs implementation
            candidate_Fmat = Fmat - rho * projected_grad;
            candidate_Fmat = (candidate_Fmat + candidate_Fmat') / 2;
            [candidate_loss, candidate_grad] = LossAndGradient(Y, As, candidate_Fmat, lambda, LR);
            
            if candidate_loss < loss || fixedStepSize
                step_size_OK = true;
                Fmat = candidate_Fmat;
                loss = candidate_loss;
                grad = candidate_grad;
                obj_seq(iter+1) = loss;
                rho_seq(iter) = rho;
                Fmat_list{iter} = Fmat;
                EY_list{iter} = AsByF(As, Fmat); % Needs implementation
                break;
            else
                rho = beta * rho;
            end
        end
        
        if rho <= 1e-8
            fprintf('step size too small, no update\n');
            break;
        end
    end
    
    % Trim the lists to the actual number of iterations
    Fmat_list = Fmat_list(~cellfun('isempty', Fmat_list));
    EY_list = EY_list(~cellfun('isempty', EY_list));
    obj_seq = obj_seq(~isinf(obj_seq));
    rho_seq = rho_seq(rho_seq ~= 0);
end

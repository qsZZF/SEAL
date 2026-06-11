function dst_val = calculate_emd_sinkhorn(dist1, dist2, cost_matrix, lambda, max_iter, tol)
%CALCULATE_EMD_SINKHORN Calculate EMD or Sinkhorn-regularized transport cost.
%
% input:
%   dist1 (Nx1 double)
%   dist2 (Nx1 double)
%   cost_matrix (NxN double): transport cost matrix
%   lambda: double, entropic regularization parameter. Smaller values are
%                       closer to exact EMD but less stable.
%   max_iter (int): Maximum number of Sinkhorn iterations. Default is 500.
%   tol (double): Convergence tolerance. Default is 1e-5.
% output:
%   dst_val (double): exact EMD for small problems when linprog is available;
%                     otherwise debiased Sinkhorn transport divergence.

    dist1 = dist1(:);
    dist2 = dist2(:);
    N = length(dist1);

    if length(dist2) ~= N
        error('dist1 and dist2 must have the same length.');
    end
    if ~isequal(size(cost_matrix), [N, N])
        error('cost_matrix must be an NxN matrix matching the distributions.');
    end
    if any(dist1 < 0) || any(dist2 < 0)
        error('Input distributions must be nonnegative.');
    end

    mass1 = sum(dist1);
    mass2 = sum(dist2);
    if mass1 <= eps || mass2 <= eps
        error('Input distributions must have positive total mass.');
    end
    if abs(mass1 - 1) > 1e-6 || abs(mass2 - 1) > 1e-6
        warning('calculate_emd_sinkhorn:NormalizeDistribution', ...
            'The sum of a distribution is not equal to 1; normalizing inputs.');
        dist1 = dist1 / mass1;
        dist2 = dist2 / mass2;
    end

    if nargin < 4 || isempty(lambda)
       lambda = 1;
    end
    if nargin < 5 || isempty(max_iter)
       max_iter = 500;
    end
    if nargin < 6 || isempty(tol)
       tol = 1e-5;
    end
    if lambda <= 0
        error('lambda must be positive.');
    end

    if N <= 200 && exist('linprog', 'file') == 2
        try
            dst_val = emd_lin(dist1, dist2, cost_matrix);
            return;
        catch ME
            warning('calculate_emd_sinkhorn:LinearProgramFailed', ...
                'Exact EMD via linprog failed (%s). Falling back to Sinkhorn.', ME.message);
        end
    end

    % Use debiased Sinkhorn divergence so identical distributions map to 0.
    ot12 = sinkhorn_ot_cost(dist1, dist2, cost_matrix, lambda, max_iter, tol);
    ot11 = sinkhorn_ot_cost(dist1, dist1, cost_matrix, lambda, max_iter, tol);
    ot22 = sinkhorn_ot_cost(dist2, dist2, cost_matrix, lambda, max_iter, tol);
    dst_val = max(ot12 - 0.5 * (ot11 + ot22), 0);
end

function dst_val = emd_lin(dist1, dist2, cost_matrix)
% Solve exact EMD using linear programming.
    N = length(dist1);
    c = cost_matrix(:);

    Aeq = zeros(2 * N, N * N);
    for i = 1:N
        % sum(F(i,:)) = dist1(i)
        indices = (0:N-1) * N + i;
        Aeq(i, indices) = 1;

        % sum(F(:,i)) = dist2(i)
        indices = (i-1) * N + (1:N);
        Aeq(N + i, indices) = 1;
    end

    beq = [dist1; dist2];
    lb = zeros(N * N, 1);
    options = optimoptions('linprog', 'Display', 'none');

    [~, dst_val, exitflag] = linprog(c, [], [], Aeq, beq, lb, [], options);
    if exitflag <= 0
        error('linprog did not converge. Exit flag: %d', exitflag);
    end
end

function ot_cost = sinkhorn_ot_cost(dist1, dist2, cost_matrix, lambda, max_iter, tol)
% Compute the entropic OT cost used by the debiased Sinkhorn divergence.
    N = length(dist1);
    K = exp(-cost_matrix / lambda);
    v = ones(N, 1);
    u = ones(N, 1);

    for i = 1:max_iter
        u_new = dist1 ./ max(K * v, realmin);
        u_new(~isfinite(u_new)) = 0;

        v_new = dist2 ./ max(K' * u_new, realmin);
        v_new(~isfinite(v_new)) = 0;

        if norm(v_new - v, 1) < tol && norm(u_new - u, 1) < tol
            u = u_new;
            v = v_new;
            break;
        end

        u = u_new;
        v = v_new;
    end

    P = bsxfun(@times, bsxfun(@times, K, u), v');
    ot_cost = sum(P(:) .* cost_matrix(:));
end

function [Source_TS, Param_TS] = seal_TS_Champagne(Data, L, VertConn, varargin)
%SEAL_TS_CHAMPAGNE Computes Two-Stage Champagne for extended source imaging.
%   [Source, Param] = seal_TS_Champagne(Data, L, VertConn, 'ParameterName', ParameterValue, ...)
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints measured data.
%       L (double matrix): Nchannels x Nsources lead field.
%       VertConn (logical matrix): Nlocations x Nlocations adjacency matrix.
%
%   Optional Name-Value Pair Inputs:
%       'NumOrientations' (integer): Orientations per location. Default: 1.
%       'PeakThreshold' (double): Rel. threshold to identify Stage 1 peaks. Default: 0.01.
%       'ExpansionHops_k' (integer): Hops to define the active neighborhood. Default: 5.
%       'BasisOrders_p' (integer): Spatial smoothing basis orders. Default: 6.
%       'ChampagneOpts' (cell): Cell array of name-value pairs for seal_Champagne.

    p = inputParser;
    p.CaseSensitive = false;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'VertConn', @(x) isnumeric(x) || islogical(x) || issparse(x));
    
    addParameter(p, 'NumOrientations', 1, @isscalar);
    addParameter(p, 'PeakThreshold', 0.01, @isscalar);
    addParameter(p, 'ExpansionHops_k', 5, @isscalar);
    addParameter(p, 'BasisOrders_p', 6, @isscalar);
    addParameter(p, 'ChampagneOpts', {}, @iscell);

    parse(p, Data, L, VertConn, varargin{:});
    opts = p.Results;

    nd = opts.NumOrientations;
    delta = opts.PeakThreshold;
    k_hops = opts.ExpansionHops_k;
    p_orders = opts.BasisOrders_p;
    
    % Inject optimal default Champagne settings if user provided none
    champ_opts = opts.ChampagneOpts;
    if isempty(champ_opts)
        champ_opts = { ...
            'MaxIterations', 100, ...
            'Tolerance', 1e-4, ...
            'PruneThreshold', 1e-8, ...
            'UpdateNoise', true ...
        };
    end

    nSource = size(L, 2);
    nLoc = size(VertConn, 1);

    %% 2. Stage 1: Focal SBL Reconstruction
    fprintf('\n--- TS-Champagne: Stage 1 (Peak Identification, nd=%d) ---\n', nd);
    
    % Pass orientations to Champagne so it groups sparsity correctly
    [S0, Param0] = seal_Champagne(Data, L, 'NumOrientations', nd, champ_opts{:});
    
    % Calculate location variance to find active anatomical peaks
    QS0 = sum(S0.^2, 2); 
    if nd > 1
        QS0_loc = sum(reshape(QS0, nd, []), 1)';
    else
        QS0_loc = QS0;
    end
    
    peak_loc_indices = find(QS0_loc > delta * max(QS0_loc));
    
    if isempty(peak_loc_indices)
        warning('seal_TS_Champagne:NoPeaks', 'Stage 1 found no peaks above threshold.');
        Source_TS = zeros(nSource, size(Data,2));
        Param_TS = struct();
        return;
    end
    fprintf('Identified %d anatomical anchor peaks.\n', length(peak_loc_indices));

    %% 3. Build Topological Neighborhood Dictionary (W2)
    fprintf('Building spatial basis dictionary...\n');
    
    % Build boolean distance matrices up to max needed order on LOCATION graph
    Q = spones(VertConn + speye(nLoc)); 
    max_order = max(k_hops + 1, p_orders + 1);
    QZ = cell(max_order, 1);
    QQ = speye(nLoc);
    
    for i = 1:max_order
        QZ{i} = double(QQ > 0); 
        if i < max_order, QQ = QQ * Q; end
    end
    
    % Find all location vertices within k_hops of the anchor peaks
    Q_neighborhood = QZ{k_hops + 1}; 
    indz_loc = [];
    for i = 1:length(peak_loc_indices)
        local_neighbors = find(Q_neighborhood(peak_loc_indices(i), :) > 0);
        indz_loc = [indz_loc, local_neighbors];
    end
    indz_loc = unique(indz_loc); 
    
    % Construct the 1D Spatial Dictionary W2_loc
    num_basis_loc = (p_orders + 1) * length(indz_loc);
    W2_loc = sparse(nLoc, num_basis_loc);
    
    col_idx = 1;
    for i = 1:(p_orders + 1)
        for j = 1:length(indz_loc)
            W2_loc(:, col_idx) = QZ{i}(:, indz_loc(j));
            col_idx = col_idx + 1;
        end
    end
    
    % Expand Spatial Dictionary to 3D Orientation Space
    if nd > 1
        W2 = kron(W2_loc, eye(nd));
    else
        W2 = W2_loc;
    end
    fprintf('Generated spatial dictionary with %d effective basis functions.\n', size(W2, 2));

    %% 4. Stage 2: Extended SBL Reconstruction
    fprintf('\n--- TS-Champagne: Stage 2 (Extended Reconstruction) ---\n');
    
    % Project Leadfield into the localized spatial dictionary
    L_dict = L * W2;
    
    % Run Champagne on the Dictionary space. (Dictionary maintains 'nd' structure)
    [S01, Param1] = seal_Champagne(Data, L_dict, 'NumOrientations', nd, champ_opts{:});
    
    % Project coefficients back to full 3D cortical space
    Source_TS = W2 * S01;

    %% 5. Output Assignment
    Param_TS.Stage1 = Param0;
    Param_TS.Stage2 = Param1;
    Param_TS.SpatialDictionary_W2 = W2;
    Param_TS.AnchorIndices = peak_loc_indices;
    Param_TS.ExpandedNeighborhoodIndices = indz_loc;
    
    fprintf('TS-Champagne Completed.\n\n');
end
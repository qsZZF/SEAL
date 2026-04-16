function [Source_MCMV, Param_MCMV] = seal_MCMV(Data, L, varargin)
%SEAL_MCMV Computes the Multi-Core Minimum Variance (MCMV) Beamformer.
%   [Source, Param] = seal_MCMV(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Generates a set of MCMV (multi-source) scalar beamformer weights. 
%   Solves the correlated-source problem by placing null-constraints on 
%   cross-talk between the specified source locations.
%
%   Based on Moiseev et al., NeuroImage, 2011.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data.
%                             Can be empty [] if DataCovariance is provided.
%       L (double matrix): Nchannels x Nsources matrix (lead field).
%                          NOTE: For MCMV, L should typically represent a 
%                          pre-selected network of interest, not the whole brain.
%
%   Optional Name-Value Pair Inputs:
%       'DataCovariance' (double matrix): Nchannels x Nchannels data covariance (C_y).
%       'NoiseCovariance' (double matrix): Nchannels x Nchannels noise covariance (C_noise).
%                                          Default: eye(Nchannels).
%       'EvokedCovariance' (double matrix): C_avg (2nd moments of epoch-averaged signals).
%                                           If provided, optimizes SNR for evoked responses.
%       'RegularizationParameter' (double): Tikhonov regularization for covariance inversion.
%                                           Default: 1e-3.
%       'NumOrientations' (integer): Number of orientations per source. Default: 1.
%
%   Outputs:
%       Source_MCMV (double matrix): Nsources x Ntimepoints estimated source activity.
%       Param_MCMV (struct):
%           .InverseOperator (double matrix): Nsources x Nchannels spatial filter.
%           .OptimalOrientations (double matrix): nd x Nlocations matrix.
%           .SourceOrder (vector): The iterative ranking of source strengths.
%           .BeamSNR (vector): The joint power SNR values.

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;
    
    defaultRegParam = 1e-3;
    defaultNumOrientations = 1;

    addRequired(p, 'Data', @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'DataCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'EvokedCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'RegularizationParameter', defaultRegParam, @isscalar);
    addParameter(p, 'NumOrientations', defaultNumOrientations, @isscalar);

    parse(p, Data, L, varargin{:});

    Data_orig = p.Results.Data;
    L_orig = p.Results.L;
    C_y = p.Results.DataCovariance;
    C_noise = p.Results.NoiseCovariance;
    C_avg = p.Results.EvokedCovariance;
    lambda = p.Results.RegularizationParameter;
    nd = p.Results.NumOrientations;

    [Nchannels, Ntimepoints] = size(Data_orig);
    [Nchannels_L, Nsources_total] = size(L_orig);
    
    if isempty(Data_orig), Nchannels = Nchannels_L; end
    if isempty(C_noise), C_noise = eye(Nchannels); end
    
    Nsources_locations = Nsources_total / nd;

    %% 2. Covariance Setup
    if isempty(C_y)
        if isempty(Data_orig), error('Must provide Data or DataCovariance.'); end
        Data_mean_sub = bsxfun(@minus, Data_orig, mean(Data_orig, 2));
        C_y = (Data_mean_sub * Data_mean_sub') / (Ntimepoints - 1);
    end
    
    % Regularize and invert Data Covariance
    C_y = (C_y + C_y') / 2;
    C_reg = C_y + lambda * (trace(C_y) / Nchannels) * eye(Nchannels);
    iR = invSPD(C_reg);

    % Scale default noise covariance if identity was passed
    if isequal(C_noise, eye(Nchannels))
        C_noise = Nchannels * eye(Nchannels) / trace(iR); 
    end

    Is_Evoked = ~isempty(C_avg);
    if Is_Evoked
        iR_working = iR * C_avg * iR;
    else
        iR_working = iR;
    end
    
    iR_N_iR = iR * C_noise * iR;

    %% 3. Independent Initialization (Single Source SNR)
    lstU3D = zeros(Nsources_locations, nd);
    SNR = zeros(Nsources_locations, 1);

    for i_loc = 1:Nsources_locations
        curIdx = (i_loc-1)*nd + 1 : i_loc*nd;
        h = L_orig(:, curIdx); 
        
        T = h' * iR_N_iR * h;
        S = h' * iR_working * h;
        
        if nd > 1
            [V, Ev] = eig(S, T);
            [~, idx] = sort(diag(Ev), 'descend');
            u = V(:, idx(1));
            u = real(u / norm(u));
        else
            u = 1;
        end
        
        lstU3D(i_loc, :) = u';
        hu = h * u;
        
        iTu = invSPD(hu' * iR_N_iR * hu);
        Su = hu' * iR_working * hu;
        SNR(i_loc) = trace(Su * iTu);
    end

    %% 4. Iterative Nulling (Moiseev 2011)
    order = 1;
    done = [];
    
    for iIter = 2:Nsources_locations
        nRef = iIter - 1;
        [~, idxMax] = sort(SNR(:), 'descend');
        order = idxMax(1:nRef);
        
        Hur = zeros(Nchannels, nRef);
        for o = 1:numel(order)
            curIdx = (order(o)-1)*nd + 1 : order(o)*nd;
            Hur(:, o) = L_orig(:, curIdx) * lstU3D(order(o), :)';
        end
        
        iR_N_iR_Hur = iR_N_iR * Hur;
        iR_Hur = iR_working * Hur;
        
        Tr = Hur' * iR_N_iR_Hur;
        Sr = Hur' * iR_Hur;
        iTr = invSPD(Tr);
        iTr_Sr_iTr = iTr * Sr * iTr;
        
        for i_loc = 1:Nsources_locations
            if any(i_loc == order) && ~any(done == i_loc)
                done = cat(1, done, i_loc);
                SNR(i_loc) = 999.9 - iIter/10; % Keep order marker
            elseif ~any(i_loc == order)
                curIdx = (i_loc-1)*nd + 1 : i_loc*nd;
                h = L_orig(:, curIdx);
                
                T = h' * iR_N_iR * h;
                S = h' * iR_working * h;
                
                Tsr = h' * iR_N_iR_Hur;
                Ssr = h' * iR_Hur;
                
                D = Tsr * iTr_Sr_iTr * Tsr' - Tsr * iTr * Ssr' - Ssr * iTr * Tsr' + S;
                F = T - Tsr * iTr * Tsr';
                
                if nd > 1
                    [V, Ev] = eig(D, F);
                    [~, idx] = sort(diag(Ev), 'descend');
                    u = V(:, idx(1));
                    u = real(u / norm(u));
                else
                    u = 1;
                end
                
                lstU3D(i_loc, :) = u';
                hu = h * u;
                Hu = horzcat(Hur, hu);
                
                iTu = invSPD(Hu' * iR_N_iR * Hu);
                Su = Hu' * iR_working * Hu;
                SNR(i_loc) = trace(iTu * Su);
            end
        end
    end

    %% 5. Compute Joint Beamformer Weights
    notRef = setdiff(1:Nsources_locations, order);
    final_order = [order(:); notRef(:)];
    
    Ht = zeros(Nsources_locations, Nchannels);
    for i = 1:Nsources_locations
        loc = final_order(i);
        curIdx = (loc-1)*nd + 1 : loc*nd;
        Ht(i, :) = (lstU3D(loc, :) * L_orig(:, curIdx)');
    end
    
    % The Joint MCMV Equation: W = (H^T * R^-1 * H)^-1 * H^T * R^-1
    wT_ordered = invSPD(Ht * iR * Ht') * Ht * iR;
    
    % Re-map weights back to original source indexing
    W_final = zeros(Nsources_locations, Nchannels);
    W_final(final_order, :) = wT_ordered;

    %% 6. Output Assignment
    Param_MCMV.InverseOperator = W_final;
    Param_MCMV.OptimalOrientations = lstU3D';
    Param_MCMV.SourceOrder = final_order';
    
    if ~Is_Evoked, SNR = SNR - Nsources_locations; end
    Param_MCMV.BeamSNR = SNR(notRef);

    if ~isempty(Data_orig)
        Source_MCMV = W_final * Data_orig;
    else
        Source_MCMV = [];
    end
end

%% Helper Function: Robust SPD Inversion
function Am1 = invSPD(A)
    try
        R = chol((A+A')/2); 
        Am1 = R \ (R' \ eye(size(A))); 
    catch
        Am1 = pinv(A);
    end
end
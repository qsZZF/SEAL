function [Source_sLORETA, Param_sLORETA] = seal_sLORETA(Data, L, varargin)
%SEAL_SLORETA Computes the standardized LORETA (sLORETA) estimate.
%   [Source_sLORETA, Param_sLORETA] = seal_sLORETA(Data, L, 'ParameterName', ParameterValue, ...)
%
%   sLORETA normalizes an MNE solution (typically with identity source covariance)
%   by a factor related to the resolution matrix M*L.
%   This implementation uses sqrt(diag(M_w * L_w')) for nd=1, and block-wise
%   sqrtm(pinv(R_block)) for nd>1.
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured EEG/MEG data.
%       L (double matrix): Nchannels x Nsources matrix (original lead field).
%
%   Optional Name-Value Pair Inputs:
%       'RegularizationParameter' (double or char/string): For the base MNE.
%                                           Default: 1/3, matching Brainstorm's
%                                           fixed SNR=3 convention for the
%                                           scaled SEAL MNE parameter.
%       'NoiseCovariance' (double matrix): Noise covariance (C). Default: eye(Nchannels).
%       'NumOrientations' (integer): Number of orientations per source location. Default: 1.
%       'SourceCovariance' (double matrix or char/string): Optional base MNE
%                                           source covariance. Default: identity.
%                                           Brainstorm disables depth weighting
%                                           for sLORETA, so the GUI does not pass
%                                           depth_weighted here.
%       % Other parameters accepted by seal_MNE can be passed.
%
%   Outputs:
%       Source_sLORETA (double matrix): Nsources x Ntimepoints sLORETA estimates.
%       Param_sLORETA (struct): Structure containing parameters and results:
%           .InverseOperator (double matrix): Operator mapping original Data to sLORETA Source.
%           .BaseMNE_Solution (double matrix): The MNE solution before normalization.
%           .BaseMNE_Parameters (struct): The 'Param' struct returned by seal_MNE.
%           .NormalizationFactors_or_Matrices (variant): Scalar factors (nd=1) or cell array of matrices (nd>1).
%           .OptionsPassed (struct): Parsed input options specific to sLORETA call.
%
%   See also: SEAL_MNE, SEAL_DSPM, RUNESIANALYSIS.

%   Author: FengZhao

    %% Input Parsing for sLORETA specific and common MNE parameters
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true; % Pass remaining to seal_MNE

    defaultRegParam_sLORETA = 1/3;
    if size(L,1) > 0
        defaultNoiseCov_sLORETA = eye(size(L, 1));
    else
        defaultNoiseCov_sLORETA = [];
    end
    defaultNumOrientations_sLORETA = 1;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

    addParameter(p, 'RegularizationParameter', defaultRegParam_sLORETA);
    addParameter(p, 'NoiseCovariance', defaultNoiseCov_sLORETA, @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
    addParameter(p, 'NumOrientations', defaultNumOrientations_sLORETA, @(x) isnumeric(x) && isscalar(x) && x>=1);
    addParameter(p, 'SourceCovariance', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x));

    try
        parse(p, Data, L, varargin{:});
    catch ME
        disp('Error parsing inputs for seal_sLORETA:');
        rethrow(ME);
    end

    Param_sLORETA.OptionsPassed = p.Results;

    current_NoiseCov = p.Results.NoiseCovariance;
    if isempty(current_NoiseCov)
        current_NoiseCov = eye(size(L,1));
        Param_sLORETA.OptionsPassed.NoiseCovariance = current_NoiseCov;
    end
    
    Nsources_total = size(L, 2);
    sourceCovariance = p.Results.SourceCovariance;
    if isempty(sourceCovariance)
        sourceCovariance = eye(Nsources_total);
    end

    mneOpts = {'SourceCovariance', sourceCovariance, ...
               'RegularizationParameter', p.Results.RegularizationParameter, ...
               'NoiseCovariance', current_NoiseCov, ...
               'NumOrientations', p.Results.NumOrientations};
               
    unmatchedArgs = namedargs2cell(p.Unmatched);
    mneOpts = [mneOpts, unmatchedArgs];

    %% 1. Compute Base MNE Solution
    [S_mne, Param_mne] = seal_MNE(Data, L, mneOpts{:});
    Param_sLORETA.BaseMNE_Solution = S_mne;
    Param_sLORETA.BaseMNE_Parameters = Param_mne;

    %% 2. Compute sLORETA Normalization and Normalized Whitened Operator M_w_normalized_rows
  
    
    M_w_mne = Param_mne.InverseOperator; 


    nd = p.Results.NumOrientations;
    M_w_normalized_rows = zeros(size(M_w_mne)); % Initialize

    if nd == 1
        resolution_diag = sum(M_w_mne .* L', 2); % Equivalent to diag(M_w_mne * L_w')
        resolution_diag(resolution_diag < 0) = 0; % Avoid sqrt of negative from numerical noise
        normalizationFactors = sqrt(resolution_diag);
        normalizationFactors(normalizationFactors < eps) = eps; % Avoid division by zero
        Param_sLORETA.NormalizationFactors_or_Matrices = normalizationFactors;
        M_w_normalized_rows = bsxfun(@rdivide, M_w_mne, normalizationFactors);
    else % nd > 1, block-wise normalization
        Nsources_locations = Nsources_total / nd;
        normalizationMatricesCell = cell(Nsources_locations, 1);
        
        for i_loc = 1:Nsources_locations
            idx = (i_loc-1)*nd + 1 : i_loc*nd;
            M_w_block = M_w_mne(idx, :);      % nd x Nchannels
            L_w_block = L(:, idx);      % Nchannels x nd
            
            R_block = M_w_block * L_w_block; % nd x nd, diagonal block of resolution matrix M_w_mne * L_w'
            
            [Ur, Sr, Vr] = svd(R_block);
            sr = diag(Sr);
            if isempty(sr) || sr(1) == 0
                rnk = 0;
            else
                rnk = sum(sr > (length(sr) * eps(single(sr(1)))));
            end

            if rnk > 0
                norm_matrix_loc = Vr(:, 1:rnk) * diag(1 ./ sqrt(sr(1:rnk))) * Ur(:, 1:rnk)';
            else
                norm_matrix_loc = zeros(nd);
            end
            normalizationMatricesCell{i_loc} = norm_matrix_loc;
            
            M_w_normalized_rows(idx, :) = norm_matrix_loc * M_w_block;
        end
        Param_sLORETA.NormalizationFactors_or_Matrices = normalizationMatricesCell;
    end

    %% 3. Compute sLORETA Estimate
    % Source_sLORETA = M_w_normalized_rows * (iW * Data); % This is Data_w
    % Or, using S_mne which is M_w_mne * Data_w
    if nd == 1
        Source_sLORETA = S_mne ./ Param_sLORETA.NormalizationFactors_or_Matrices;
    else % For nd > 1, M_w_normalized_rows already incorporates the normalization transform
        Source_sLORETA = M_w_normalized_rows  * Data;
    end


    %% 4. Compute sLORETA Inverse Operator
    Param_sLORETA.InverseOperator = M_w_normalized_rows;

end

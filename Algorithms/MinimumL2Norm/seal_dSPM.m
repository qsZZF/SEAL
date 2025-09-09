function [Source_dSPM, Param_dSPM] = seal_dSPM(Data, L, varargin)
%SEAL_DSPM Computes the dynamic Statistical Parametric Mapping (dSPM).
%   [Source_dSPM, Param_dSPM] = seal_dSPM(Data, L, 'ParameterName', ParameterValue, ...)
%
%   dSPM normalizes an MNE solution (typically with identity source covariance)
%   by a factor related to the noise sensitivity of each source estimate.
%   This implementation uses row norms of the MNE inverse operator for normalization,
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured EEG/MEG data.
%       L (double matrix): Nchannels x Nsources matrix (original lead field).
%
%   Optional Name-Value Pair Inputs:
%       'RegularizationParameter' (double or char/string): regularization value for Tikhonov regularization 
%                                           lambda (the value used in seal_MNE is lambda^2*trace(L*L')/trace(Noise_Cov)).
%                                           If 'auto', attempts to estimate from data.
%                                           Default: 0.1
%       'NoiseCovariance' (double matrix): Nchannels x Nchannels noise covariance
%                                          matrix (C). Default: eye(Nchannels).
%       'NumOrientations' (integer): Number of orientations per source location (e.g., 1 for fixed, 3 for free).
%                                    If > 1, L should be Nchannels x (Nsources_locations * NumOrientations).
%                                    SourceCovariance 'auto' will produce a diagonal spatial covariance
%                                    which is then kron'd with eye(NumOrientations).
%                                    Default: 1.
%       Other parameters accepted by seal_MNE can be passed. Typically there is no
%       need to apply depth compensation 
%
%   Outputs:
%       Source_dSPM (double matrix): Nsources x Ntimepoints dSPM estimates.
%       Param_dSPM (struct): Structure containing parameters and results:
%           .InverseOperator (double matrix): Nsources x Nchannels inverse operator 
%                                             to be applied to original (unwhitened) data.
%           .BaseMNE_Solution (double matrix): The MNE solution before normalization.
%           .BaseMNE_Parameters (struct): The 'Param' struct returned by seal_MNE.
%           .NormalizationFactors (double vector): Nsources x 1 normalization factors.
%           .OptionsPassed (struct): Parsed input options specific to dSPM call.


%   See also: SEAL_MNE, SEAL_SLORETA,

%   Author: FengZhao

    %% Input Parsing for dSPM specific and common MNE parameters
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true; % Pass remaining to seal_MNE

    % Parameters for seal_MNE (can be overridden by user)
    defaultRegParam_dSPM = 0.1; % Default for the base MNE
    defaultNoiseCov_dSPM = eye(size(L, 1));
    defaultNumOrientations_dSPM = 1;

    addRequired(p, 'Data', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));

    % Explicitly define parameters that dSPM might want to control for seal_MNE,
    % or dSPM-specific parameters if any. For now, mainly pass-through.
    addParameter(p, 'RegularizationParameter', defaultRegParam_dSPM);
    addParameter(p, 'NoiseCovariance', defaultNoiseCov_dSPM, @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'NumOrientations', defaultNumOrientations_dSPM, @(x) isnumeric(x) && isscalar(x) && x>=1);

    try
        parse(p, Data, L, varargin{:});
    catch ME
        disp('Error parsing inputs for seal_dSPM:');
        rethrow(ME);
    end

    % Store parsed options for dSPM to record what was passed to it
    Param_dSPM.OptionsPassed = p.Results;
    
    % Options to pass to seal_MNE
    % dSPM/sLORETA typically use identity source covariance for the base MNE
    Nsources_total = size(L, 2);
    mneOpts = {'SourceCovariance', eye(Nsources_total), ...
               'RegularizationParameter', p.Results.RegularizationParameter, ...
               'NoiseCovariance', p.Results.NoiseCovariance, ...
               'NumOrientations', p.Results.NumOrientations};
           
    % Append any unmatched parameters from varargin for seal_MNE
    unmatchedArgs = namedargs2cell(p.Unmatched);
    mneOpts = [mneOpts, unmatchedArgs];

    %% 1. Compute Base MNE Solution
    [S_mne, Param_mne] = seal_MNE(Data, L, mneOpts{:});
    Param_dSPM.BaseMNE_Solution = S_mne;
    Param_dSPM.BaseMNE_Parameters = Param_mne;

    M_for_norm = Param_mne.InverseOperator; % This is M_w without whitening
    %% 2. Compute dSPM Normalization Factors



    nd = p.Results.NumOrientations;
    if nd == 1
        row_wise_sum_sq = sum(M_for_norm.^2, 2);
    else 
        row_wise_sum_sq_temp = sum(M_for_norm.^2, 2); % Nsources_total x 1
        Nsources_locations = Nsources_total / nd;
        
        % Reshape into NdxNloc, sum along Ndim (dim 1), then repmat
        sum_sq_norms_at_loc = sum(reshape(row_wise_sum_sq_temp, nd, Nsources_locations), 1); % 1 x Nloc
        
        % Repeat this sum for each orientation at that location
        row_wise_sum_sq = repelem(sum_sq_norms_at_loc, nd)'; % Nsources_total x 1
    end
    
    normalizationFactors = sqrt(row_wise_sum_sq);
    normalizationFactors(normalizationFactors < eps) = eps; % Avoid division by zero

    Param_dSPM.NormalizationFactors = normalizationFactors;
    Param_dSPM.InverseOperator = bsxfun(@rdivide, M_for_norm, normalizationFactors);

    %% 3. Compute dSPM Estimate
    Source_dSPM = S_mne ./ normalizationFactors; % Broadcasting division

end

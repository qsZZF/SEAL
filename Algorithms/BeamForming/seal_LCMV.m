function [Source_LCMV, Param_LCMV] = seal_LCMV(Data, L, varargin)
%SEAL_LCMV Computes Linearly Constrained Minimum Variance (LCMV) Beamformer.
%   [Source, Param] = seal_LCMV(Data, L, 'ParameterName', ParameterValue, ...)
%
%   Inputs:
%       Data (double matrix): Nchannels x Ntimepoints matrix of measured data.
%                             Can be empty [] if DataCovariance is provided.
%       L (double matrix): Nchannels x Nsources matrix (lead field).
%
%   Optional Name-Value Pair Inputs:
%       'Type' (char/string): Beamformer variant:
%                             'ug'  - Unit-Gain (Standard LCMV)
%                             'ua'  - Unit-Array (Leadfield normalized)
%                             'ung' - Unit-Noise-Gain (Pseudo-Z / NAI)
%                             Default: 'ug'
%       'DataCovariance' (double matrix): Nchannels x Nchannels data covariance (C_y).
%       'NoiseCovariance' (double matrix): Nchannels x Nchannels noise covariance (C_noise).
%                                          If provided, data, leadfield, and data covariance
%                                          are spatially whitened prior to filter computation.
%                                          Default: eye(Nchannels).
%       'RegularizationParameter' (double): Tikhonov regularization (shrinkage) applied 
%                                           to the data covariance matrix. 
%                                           Default: 1e-3
%       'NumOrientations' (integer): Number of orientations per source. Default: 1.
%       'ScalarBeamformer' (boolean): If true and NumOrientations > 1, projects the 
%                                     vector beamformer onto the optimal orientation.
%
%   Outputs:
%       Source_LCMV (double matrix): Nsources x Ntimepoints estimated source activity.
%       Param_LCMV (struct):
%           .InverseOperator (double matrix): The Nsources x Nchannels spatial filter 
%                                             (maps original unwhitened data to source).
%           .WhiteningMatrix (double matrix): The iW applied to the system.
%           .DataCovariance (double matrix): The regularized, whitened data covariance used.
%           .OptimalOrientations (double matrix): nd x Nlocations matrix of optimal dipoles.

    %% 1. Input Parsing
    p = inputParser;
    p.CaseSensitive = false;
    
    defaultType = 'ug';
    defaultRegParam = 1e-3;
    defaultNumOrientations = 1;
    defaultScalar = false;

    addRequired(p, 'Data', @(x) (isnumeric(x) && ismatrix(x)) || isempty(x));
    addRequired(p, 'L', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'Type', defaultType, @(x) ischar(x) || isstring(x));
    addParameter(p, 'DataCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'NoiseCovariance', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
    addParameter(p, 'RegularizationParameter', defaultRegParam, @isscalar);
    addParameter(p, 'NumOrientations', defaultNumOrientations, @isscalar);
    addParameter(p, 'ScalarBeamformer', defaultScalar, @islogical);

    parse(p, Data, L, varargin{:});

    Param_LCMV.OptionsPassed = p.Results;
    Data_orig = p.Results.Data;
    L_orig = p.Results.L;
    bf_type = lower(p.Results.Type);
    C_y = p.Results.DataCovariance;
    C_noise = p.Results.NoiseCovariance;
    lambda = p.Results.RegularizationParameter;
    nd = p.Results.NumOrientations;
    use_scalar = p.Results.ScalarBeamformer;

    [Nchannels_L, Nsources_total] = size(L_orig);
    Nchannels = Nchannels_L;
    
    if ~isempty(Data_orig) && size(Data_orig, 1) ~= Nchannels
        error('seal_LCMV:DimensionMismatch', 'Data and L must have same number of channels.');
    end
    if isempty(C_noise), C_noise = eye(Nchannels); end
    
    Nsources_locations = Nsources_total / nd;

    %% 2. Noise Covariance Whitening (SVD Based)
    iW = eye(Nchannels); 
    if ~isequal(C_noise, eye(Nchannels))
        try  
            [U_c, S_c, ~] = svd((C_noise + C_noise')/2);
            s_c = diag(S_c);
            quan_c = mean(s_c)*1e-6;
            s_c(s_c < quan_c) = quan_c; % Robust thresholding
            iW = diag(1./sqrt(s_c)) * U_c';
        catch ME_whitening
            warning('seal_LCMV:WhiteningFailed', 'Whitening failed. Using identity. Error: %s', ME_whitening.message);
            iW = eye(Nchannels);
        end
    end
    Param_LCMV.WhiteningMatrix = iW;

    % Project Leadfield and Data into Whitened Space
    L_w = iW * L_orig;
    if ~isempty(Data_orig)
        Data_w = iW * Data_orig;
    else
        Data_w = [];
    end

    %% 3. Data Covariance Estimation & Whitening
    if isempty(C_y)
        if isempty(Data_w)
            error('seal_LCMV:NoCovariance', 'Must provide Data or DataCovariance.');
        end
        % Compute directly from whitened data
        Ntimepoints = size(Data_w, 2);
        Data_mean_sub = bsxfun(@minus, Data_w, mean(Data_w, 2));
        C_y_w = (Data_mean_sub * Data_mean_sub') / (Ntimepoints - 1);
    else
        % Whiten the user-provided data covariance: C_yw = iW * Cy * iW^T
        C_y_w = iW * C_y * iW';
    end
    
    % Ensure symmetry and apply Tikhonov regularization (Shrinkage)
    C_y_w = (C_y_w + C_y_w') / 2;
    C_reg = C_y_w + lambda * (trace(C_y_w) / Nchannels) * eye(Nchannels);
    Param_LCMV.DataCovariance = C_reg;

    % Inverse Covariance
    if rcond(C_reg) < 1e-12
        invC = pinv(C_reg);
    else
        invC = inv(C_reg);
    end

    %% 4. Compute Beamformer Weights (in Whitened Space)
    if use_scalar && nd > 1
        W_w_final = zeros(Nsources_locations, Nchannels);
    else
        W_w_final = zeros(Nsources_total, Nchannels);
    end
    optO = zeros(nd, Nsources_locations);

    for i_loc = 1:Nsources_locations
        curIdx = (i_loc-1)*nd + 1 : i_loc*nd;
        L_k = L_w(:, curIdx); % Use whitened leadfield
        
        % Optimal Orientation Calculation
        if nd > 1
            switch bf_type
                case {'ug', 'ua'}
                    [v, d] = eig(L_k' * invC * L_k);
                    [~, sortIdx] = sort(diag(d), 'ascend'); 
                case 'ung'
                    [v, d] = eig(L_k' * invC * L_k, L_k' * (invC * invC) * L_k);
                    [~, sortIdx] = sort(diag(d), 'ascend'); 
            end
            optO(:, i_loc) = v(:, sortIdx(1));
        else
            optO(:, i_loc) = 1;
        end

        if use_scalar && nd > 1
            L_k = L_k * optO(:, i_loc); 
        end
        
        switch bf_type
            case 'ug' % Unit Gain
                Term_inv = pinv(L_k' * invC * L_k);
                ft = Term_inv * (L_k' * invC);
                
            case 'ua' % Unit Array
                norm_factors = sqrt(sum(L_k.^2, 1));
                fnt = bsxfun(@rdivide, L_k, norm_factors);
                Term_inv = pinv(fnt' * invC * fnt);
                ft = Term_inv * (fnt' * invC);
                
            case 'ung' % Unit Noise Gain (NAI)
                Term_inv = pinv(L_k' * invC * L_k);
                w_ug = Term_inv * (L_k' * invC); 
                % Because space is whitened, noise is I, so gamma is just sum(w_ug.^2)
                gamma_diag = sum(w_ug.^2, 2); 
                ft = bsxfun(@times, w_ug, 1./sqrt(gamma_diag));
                
            otherwise
                error('seal_LCMV:UnknownType', 'Unknown LCMV Type.');
        end
        
        if use_scalar && nd > 1
            W_w_final(i_loc, :) = ft;
        else
            W_w_final(curIdx, :) = ft;
        end
    end

    %% 5. Un-whiten and Apply Spatial Filter
    % M_final = W_whitened * iW
    Param_LCMV.InverseOperator = W_w_final * iW;
    Param_LCMV.OptimalOrientations = optO;

    if ~isempty(Data_orig)
        Source_LCMV = Param_LCMV.InverseOperator * Data_orig;
    else
        Source_LCMV = [];
    end

end
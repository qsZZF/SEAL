function reg = seal_getAlgorithmRegistry()
% 每个算法 = 一份配置；新增算法只改这里
% parent: 'minL2' | 'spatio' | 'sparse'
% params: 每个参数 = {字段名, 显示标签, 控件类型, 默认值, [选项]}

reg = struct();

reg.MNE = struct( ...
    'parent','minL2', ...
    'title','MNE Settings', ...
    'params',{{
    {'RegularizationParameter','RegParam','numeric',1/3}
    {'sourceCovariance','SourceCov','edit','brainstorm_depth'}
    }});

reg.dSPM = struct('parent','minL2','title','dSPM Settings','params',{{
    {'RegularizationParameter','RegParam','numeric',0.1}
    }});

reg.sLORETA = struct('parent','minL2','title','sLORETA Settings','params',{{
    {'RegularizationParameter','RegParam','numeric',1/3}
    }});

reg.eLORETA = struct('parent','minL2','title','eLORETA Settings','params',{{ ...
    {'RegularizationParameter','RegParam','numeric',0.1}...
    {'sourceCovariance','SourceCov','edit','depth_weighted'} ...
    {'MaxIterations','MaxIterations','numeric',100} ...
    {'Tolerance','Tolerance','numeric', 1e-6} ...
    }});

reg.LORETA  = struct('parent','minL2','title','LORETA Settings','params',{{ ...
    {'RegularizationParameter','RegParam','numeric',0.1} ...
    {'DepthWeightingExponent','DepthWeightingExponent','numeric',1} ...
    {'LaplacianRegFactor','LaplacianRegFactor','numeric',1e-4} ...
    {'sourceCovariance','SourceCov','edit','depth_weighted'} ...
    }});

reg.LAURA   = struct('parent','minL2','title','LAURA Settings','params',{{ ...
    {'RegularizationParameter','RegParam','numeric',0.1} ...
    {'DepthWeightingExponent','DepthWeightingExponent','numeric',1} ...
    {'Exponent_e_LAURA','Exponent_e_LAURA','numeric',2} ...
    {'MaxNeighborsN_LAURA','MaxNeighborsN_LAURA','numeric',9}
    }});

reg.LOGURA   = struct('parent','minL2','title','LOGURA Settings','params',{{ ...
    {'RegularizationParameter','RegParam','numeric',0.1} ...
    {'DepthWeightingExponent','DepthWeightingExponent','numeric',1} ...
    {'Exponent_e_LOGURA','Exponent_e_LOGURA','numeric',2} ...
    }});



reg.BlockChampagne = struct( ...
    'parent','sparse', ...
    'title','BlockChampagne Settings', ...
    'params',{{
    {'AtlasConstraints','AtlasConstraints','numeric',0}
    {'NeighborOrder','NeiborOrder','numeric',1}
    {'LearnNoise','LearnNoise','dropdown','true',{'true','false'}}
    {'NoiseObservationMatrix','NoiseObservation','dropdown','eye',{'eye'}}
    }});

reg.Champagne = struct( ...
    'parent','sparse', ...
    'title','Champagne Settings', ...
    'params',{{
    {'PruneThreshold','PruneThreshold','numeric',1e-6}
    {'UpdateNoise','UpdateNoise','dropdown','true',{'true','false'}}
    {'MaxIterations','MaxIterations','numeric',100}
    {'Tolerance','Tolerance','numeric', 1e-8}
    }});

reg.TS_Champagne = struct( ...
    'parent','sparse', ...
    'title','TS_Champagne Settings', ...
    'params',{{ ...
     {'PeakThreshold','PeakThreshold','numeric',0.01}
    {'ExpansionHops_k','ExpansionHops_k','numeric',5}
    {'BasisOrders_p','BasisOrders_p','numeric',6}
    }});
%  'customBuilder','No need');

reg.SmoothChampagne = struct( ...
    'parent','sparse', ...
    'title','Smooth-Cham Settings', ...
    'params',{{ ...
    {'KernelWidth','KernelWidth (m)','numeric',2}
    {'KernelPower','KernelPower (p)','numeric',4}
    {'KernelThreshold','KernelThreshold','numeric',0.1}
    {'TileSize','TileSize','numeric',6}
    {'MaxIterations','MaxIterations','numeric',100}
    {'Tolerance','Tolerance','numeric',1e-6}
    {'UpdateNoiseCovariance','UpdateNoiseCovariance','dropdown','false',{'true','false'}}
    {'Verbose','Verbose','dropdown','false',{'true','false'}}
    }});

reg.IRES = struct( ...
    'parent','sparse', ...
    'title','IRES Settings', ...
    'params',{{ ...
    {'Alpha','Alpha','numeric',0.05}
    {'MaxIterations','ReweightIterations','numeric',4}
    {'MaxADMMIterations','MaxADMMIterations','numeric',200}
    {'Tolerance','ReweightTolerance','numeric',1e-3}
    {'ADMMTolerance','ADMMTolerance','numeric',1e-4}
    {'ProgressInterval','ProgressInterval','numeric',25}
    {'Rho','ADMMRho','numeric',100}
    {'WeightEpsilon','WeightEpsilon','numeric',1e-3}
    {'Verbose','Verbose','dropdown','true',{'true','false'}}
    }});

reg.FOCUSS = struct( ...
    'parent','sparse', ...
    'title','FOCUSS Settings', ...
    'params',{{ ...
    {'RegularizationParameter','RegularizationParameter','numeric',0.005}
    {'RegularizationMode','RegularizationMode','dropdown','fixed_initial',{'fixed_initial','adaptive_legacy'}}
    {'MaxIterations','MaxIterations','numeric',40}
    {'Tolerance','Tolerance','numeric',1e-4}
    {'PruneThreshold','PruneThreshold','numeric',1e-6}
    {'DepthWeighting','DepthWeighting','dropdown','none',{'none','leadfield_norm'}}
    {'WeightMode','WeightMode','dropdown','compound',{'compound','simple'}}
    {'Verbose','Verbose','dropdown','false',{'true','false'}}
    }});

reg.Lq_ADMM = struct( ...
    'parent','sparse', ...
    'title','Lq_ADMM Settings', ...
    'params',{{
    {'RegularizationParameter','RegParam','numeric',0.1}...
    {'DataFidelity','DataFidelity','edit','L2'}...
    {'ADMM_Rho','ADMM_Rho','numeric', 1}...
    {'q_norm','q_norm','numeric',0.5} ...
    {'MaxIterations','MaxIterations','numeric',500} ...
    {'Tolerance','Tolerance','numeric', 1e-5} ...
    }});
%   'customBuilder','Phi');

reg.uSTAR = struct( ...
    'parent','spatio', ...
    'title','uSTAR Settings', ...
    'params',{{
    {'NumTBFs','NumTBFs','numeric',4}
    {'MaxIterations','MaxIterations','numeric',50}
    {'Tolerance','Tolerance','numeric',1e-6}
    {'uSTARRunAlgorithm','RunAlgorithm','dropdown','STAR',{'STAR','sSTARTS'}}
    {'PriorType','PriorType','dropdown','Laplacian',{'Laplacian','LAURA'}}
    {'CwkRegularization','CwkRegularization','numeric',1e-10}
    % 特殊控件单独标记，引擎跳过，由专用方法补
    }});
% 'customBuilder','builduSTARExtras');   % 处理 MicrostateMaps 这种非标准控件


reg.STARTS = struct( ...
    'parent','spatio', ...
    'title','STARTS Settings', ...
    'params',{{
    {'NumTBFs','NumTBFs','numeric',4}
    {'MaxIterations','MaxIterations','numeric',100}
    {'Tolerance','Tolerance','numeric',1e-6}
    {'UpdateNoiseCovariance','UpdateNoiseCovariance(0/1/2)','numeric',1}
    {'NeighborOrder','NeighborOrder','numeric', 1}
    {'BlockCorrelation','BlockCorrelation','numeric',0.99} 
    {'PruneGammas','PruneGammas','numeric',1e-4} 
    {'PruneAlphas','PruneAlphas','numeric', 1e-4} 
    {'CwkRegularization','CwkRegularization','numeric',1e-10}
    }});
%InitialTBFs

reg.sSTARTS = struct( ...
    'parent','spatio', ...
    'title','sSTARTS Settings', ...
    'params',{{
    {'NumTBFs','NumTBFs','numeric',4}
    {'LaplacianTau','LaplacianTau','numeric',0.99}
    {'MaxIterations','MaxIterations','numeric',50}
    {'Tolerance','Tolerance','numeric',1e-6}
    {'UpdateNoiseCovariance','UpdateNoiseCovariance(0/1/2)','numeric',1}
    {'PruneGammas','PruneGammas','numeric',1e-6}
    {'PruneAlphas','PruneAlphas','numeric',0.1}
    {'CwkRegularization','CwkRegularization','numeric',1e-10}
    {'Verbose','Verbose','dropdown','true',{'true','false'}}
    }});
%InitialTBFs

reg.BESTIES = struct( ...
    'parent','spatio', ...
    'title','BESTIES Settings', ...
    'params',{{
    {'NumTBFs','NumTBFs','numeric',4}
    {'TemporalBasisMode','TemporalBasisMode','dropdown','svd',{'svd','none'}}
    {'InitialTau','InitialTau','numeric',0.98}
    {'LearnTau','LearnTau','dropdown','true',{'true','false'}}
    {'ComputeFreeEnergy','ComputeFreeEnergy','dropdown','true',{'true','false'}}
    {'NoiseObservationMatrix','NoiseObservation','dropdown','eye',{'eye'}}
    {'MaxIterations','MaxIterations','numeric',200}
    {'Tolerance','Tolerance','numeric', 1e-6}
    }});

reg.STBFSI = struct( ...
    'parent','spatio', ...
    'title','SI-STBF Settings', ...
    'params',{{
    {'NumTBFs','NumTBFs','numeric',4}
    {'MaxSpatialBases','MaxSpatialBases','numeric',120}
    {'NeighborOrder','NeighborOrder','numeric',1}
    {'DiffusionSigma','DiffusionSigma','numeric',0.6}
    {'DiffusionOrder','DiffusionOrder','numeric',8}
    {'SeedSelection','SeedSelection','dropdown','mne',{'mne','leadfield','uniform'}}
    {'MaxIterations','MaxIterations','numeric',50}
    {'Tolerance','Tolerance','numeric',1e-6}
    {'PruneAlphas','PruneAlphas','numeric',1e-4}
    {'PruneGammas','PruneGammas','numeric',1e-4}
    {'PruneSpatialBases','PruneSpatialBases','dropdown','true',{'true','false'}}
    {'CovarianceMode','CovarianceMode','dropdown','auto',{'auto','full','diagonal'}}
    {'Verbose','Verbose','dropdown','true',{'true','false'}}
    }});

reg.MxNE = struct( ...
    'parent','spatio', ...
    'title','MxNE Settings', ...
    'params',{{
    {'Alpha','Alpha','numeric',15}
    {'WeightsMin','WeightsMin','numeric',0.0}
    {'SolverTolerance','SolverTolerance','numeric',1e-4}
    {'DebiasMaxIter','DebiasMaxIter','numeric',1000}
    {'DebiasTolerance','DebiasTolerance','numeric',1e-7}
    {'SolverMaxIter','SolverMaxIter','numeric',3000}
    {'Debias','Debias','dropdown','false',{'true','false'}}
    {'PCAWhitening','PCAWhitening','dropdown','true',{'true','false'}}
    }});
%DepthBiasComp,PickOri,Weights

reg.irMxNE = struct( ...
    'parent','spatio', ...
    'title','irMxNE Settings', ...
    'params',{{
    {'Alpha','Alpha','numeric',10}
    {'NumIRMXNEIter','NumIRMXNEIter','numeric',15}
    {'WeightsMin','WeightsMin','numeric',0.0}
    {'SolverTolerance','SolverTolerance','numeric',1e-4}
    {'IRMXNETolerance','IRMXNETolerance','numeric',1e-6}
    {'DebiasMaxIter','DebiasMaxIter','numeric',1000}
    {'DebiasTolerance','DebiasTolerance','numeric',1e-7}
    {'SolverMaxIter','SolverMaxIter','numeric',3000}
    {'Debias','Debias','dropdown','false',{'true','false'}}
    {'PCAWhitening','PCAWhitening','dropdown','true',{'true','false'}}
    }});
%DepthBiasComp,RankNoiseCov,PickOri,Weights


reg.FASTIRES = struct( ...
    'parent','spatio', ...
    'title','FASTIRES Settings', ...
    'params',{{
    {'Alpha','Alpha','numeric',0.08}
    {'NumTBFs','NumTBFs','numeric',4}
    {'Epsilon','Epsilon','numeric',1e-3}
    {'MaxIterations','MaxIterations','numeric',10}
    {'MaxInnerIterationsX','MaxInnerIterationsX','numeric',20}
    {'MaxInnerIterationsY','MaxInnerIterationsY','numeric',20}
    {'MaxInnerIterationsProj','MaxInnerIterationsProj','numeric',20}
    {'Tolerance','Tolerance','numeric',1e-4}
    }});%InitialTBFs

reg.VSSI_Lp = struct( ...
    'parent','spatio', ...
    'title','VSSI-Lp Settings', ...
    'params',{{ ...
    {'RegularizationParameter','RegularizationParameter','numeric',0.05}
    {'P','LpExponent','numeric',0.6}
    {'Rho','Rho','numeric',1.0}
    {'MaxIterations','MaxIterations','numeric',500}
    {'Tolerance','Tolerance','numeric',1e-4}
    {'ProgressInterval','ProgressInterval','numeric',25}
    {'AdaptiveRho','AdaptiveRho','dropdown','true',{'true','false'}}
    {'NormalizeProblem','NormalizeProblem','dropdown','true',{'true','false'}}
    {'Verbose','Verbose','dropdown','true',{'true','false'}}
    }});

reg.VSSI_GGD = struct( ...
    'parent','spatio', ...
    'title','VSSI-GGD Settings', ...
    'params',{{ ...
    {'SparseWeight','SparseWeight','numeric',0.01}
    {'P','GGDShape','numeric',0.7}
    {'MaxOuterIterations','MaxOuterIterations','numeric',10}
    {'MaxADMMIterations','MaxADMMIterations','numeric',400}
    {'MaxIRLSIterations','MaxIRLSIterations','numeric',10}
    {'Tolerance','Tolerance','numeric',1e-3}
    {'ADMMTolerance','ADMMTolerance','numeric',1e-5}
    {'ProgressInterval','ProgressInterval','numeric',25}
    {'Rho','Rho','numeric',1.0}
    {'HutchinsonProbes','HutchinsonProbes','numeric',8}
    {'AdaptiveRho','AdaptiveRho','dropdown','true',{'true','false'}}
    {'NormalizeProblem','NormalizeProblem','dropdown','true',{'true','false'}}
    {'Verbose','Verbose','dropdown','true',{'true','false'}}
    }});

reg.STOUT = struct( ...
    'parent','spatio', ...
    'title','STOUT Settings', ...
    'params',{{ ...
    {'ReferenceSpectrum','ReferenceSpectrum','edit','motor'}
    {'ConstraintMode','ConstraintMode','dropdown','hyco',{'hyco','coco','exco'}}
    {'HybridGamma','HybridGamma','numeric',10}
    {'RegularizationParameter','RegularizationParameter','numeric',0.05}
    {'Verbose','Verbose','dropdown','false',{'true','false'}}
    }});

reg.dMAP = struct( ...
    'parent','spatio', ...
    'title','dMAP Settings', ...
    'params',{{
    {'lambda_F','lambda_F','numeric',0.95}
    {'alpha_F','alpha_F','numeric',0.51}
    {'max_iter','max_iter','numeric',15}
    {'tol','tol','numeric',1e-5}
    }});%InitialSourceCovDiag

reg.dSPN = struct( ...
    'parent','spatio', ...
    'title','dSPN Settings', ...
    'params',{{
    {'InitialSourceCovScale','InitialSourceCovScale','numeric',1e-3}
    {'max_iter','max_iter','numeric',30}
    {'tol','tol','numeric',1e-5}
    {'rho','rho','numeric',0.99}
    {'q_min','q_min','numeric',1e-10}
    {'sigma_min','sigma_min','numeric',1e-6}
    {'riccati_iter','riccati_iter','numeric',200}
    {'riccati_tol','riccati_tol','numeric',1e-5}
    {'InnovationJitter','InnovationJitter','numeric',1e-8}
    {'smooth_pass','smooth_pass','dropdown','true',{'true','false'}}
    }});


reg.debiased_GroupLasso = struct( ...
    'parent','spatio', ...
    'title','debiased_GroupLasso Settings', ...
    'params',{{
    {'RegularizationParameter','RegParam','numeric',0.1}
    {'DynamicLambda','DynamicLambda','dropdown','true',{'true','false'}}
    {'JointEstNum','JointEstNum','numeric',10}
    {'Debias','Debias','dropdown','true',{'true','false'}}
    {'MaxIterations','MaxIterations','numeric',1000}
    {'Tolerance','Tolerance','numeric',1e-8}
    {'ToleranceNorm','ToleranceNorm','numeric',1e-3}
    {'JointTolerance','JointTolerance','numeric',1e-4}
    {'VaryingRho','VaryingRho','dropdown','true',{'true','false'}}
    {'ClearNotSelect','ClearNotSelect','dropdown','true',{'true','false'}}
    {'PrewhitenNoise','PrewhitenNoise','dropdown','true',{'true','false'}}
    }});

reg.SISSY_L21 = struct( ...
    'parent','spatio', ...
    'title','SISSY_L21 Settings', ...
    'params',{{
    {'lambda','lambda','numeric',1.0}
    {'DynamicLambda','DynamicLambda','dropdown','true',{'true','false'}}
    {'JointEstNum','JointEstNum','numeric',10}
    {'alpha','alpha','numeric',0.07}
    {'rho','rho','numeric',1.0}
    {'max_iter','max_iter','numeric',200}
    {'tol','tol','numeric',1e-4}
    }});

reg.TFMxNE = struct( ...
    'parent','spatio', ...
    'title','TFMxNE Settings', ...
    'params',{{
    {'SpatialReg','SpatialReg','numeric',50}
    {'TemporalReg','TemporalReg','numeric',1}
    {'TimeStep','TimeStep','numeric',4}
    {'WindowSize','WindowSize','numeric',64}
    {'MaxIterations','MaxIterations','numeric',200}
    {'Tolerance','Tolerance','numeric',1e-4}
    {'DepthBiasExponent','DepthBiasExponent','numeric',0.8}
    {'DepthBiasLimit','DepthBiasLimit','numeric',10}
    {'LooseOrientationWeight','LooseOriWeight','numeric',1.0}
    {'Debias','Debias','dropdown','false',{'true','false'}}
    }});

reg.MCMV = struct( ...
    'parent','beamformer', ...
    'title','MCMV Settings', ...
    'params',{{
    {'RegularizationParameter','RegParam','numeric',0.05}
    {'MaxIterations','MaxIterations','numeric',30}
    {'AutoScaleNoise','AutoScaleNoise','dropdown','true',{'true','false'}}
    {'Verbose','Verbose','dropdown','false',{'true','false'}}
    }});%DataCovariance,EvokedCovariance

reg.LCMV = struct( ...
    'parent','beamformer', ...
    'title','LCMV Settings', ...
    'params',{{
       {'Type','Type','edit','ug'}
        {'RegularizationParameter','RegParam','numeric',0.05}
        {'ScalarBeamformer','ScalarBeamformer','dropdown','true',{'true','false'}}
    }});%DataCovariance


end

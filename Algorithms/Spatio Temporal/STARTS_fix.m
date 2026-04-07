function [S,par] = STARTS_fix(L,B,Phi,MM,varargin)
% This script is to solve E/MEG source imaging problem under STARTS
% framework
% Source Model: S = W*Phi, B = LS+E+epsilon
% =================================================================
% [S,par] = STARTS_fix(L,B,Phi,M,varargin)
% Input:
%   -L: lead-field matrix
%   -B: E/MEG measurements
%   -Phi: temporal basis functions (TBFs)
%   -M: spatial prior
% Output:
%   -S: reconstructed source
%   -par:
%       par.Phi
%       par.W
%       par.E
%==================================================================
% Author: Feng Zhao
% Data: 2022/5/8
% Reference:
% [1] STARTS: A Self-Adapted Spatio-Temporal Framework for Automatic E/MEG Source Imaging

%% check input
maxIter = 50;
epsilon = 1e-6;
updateR = 0;
Knb = 1; % Kth of dipoles
isextended = 1;
if nargin>4
    for inar = 1:2:length(varargin)
        Param = lower(varargin{inar});
        Value = varargin{inar+1};
        switch Param
            case 'maxiter'
                maxIter = Value;
            case 'epsilon'
                epsilon = Value;
            case 'updater'
                updateR = Value;
            case 'knb'
                Knb = Value;
            case 'isextended'
                isextended = Value;
        end
    end
end
%% Initialize
[nSensor,nSource] = size(L); % number of channels & source dipoles

nSamp = size(B,2); % snapshots
if norm(B,'fro')<sqrt(numel(B))
    B = B./norm(B,'fro')*sqrt(numel(B));
end
if isempty(Phi)
    [~,tS,tD] = svd(B);
    tS = diag(tS);
    [tS,tsIdx] = sort(tS,'descend');
    tD = tD(:,tsIdx);
    numPhi = max(1,find(cumsum(tS)./sum(tS)>0.9,1));
    numPhi = 4;
    Phi = tD(:,1:numPhi)';
end
beta = 1e-4;
C_noise = eye(nSensor); % noise covariance

if isextended
    LF = [L eye(nSensor)]; % leadfield: source & noise gain
    naSource = nSource+nSensor;
    M = blkdiag(MM^Knb,zeros(nSensor));
    M = M+eye(size(M));
else
    LF = L; % leadfield: source & noise gain
    naSource = nSource;
    M = MM^Knb;
    M = M+eye(size(M));
end
cls = sum(M>0,2);
G = sparse(naSource,sum(cls));
clusters = cell(numel(cls),1);
clusterBlk = clusters;

Rcorr = 0.99; % inter-block correlation coefficients

for it = 1:numel(cls)
    if it == 1
        clusters{it} = 1:sum(cls(1:it));
    else
        clusters{it} = sum(cls(1:it-1))+1 : sum(cls(1:it));
    end
    tidx = sort(unique([it,find(M(it,:)>0)]));
    if cls(it) == 1
        clusterBlk{it} = 1;
        G(tidx,clusters{it}) = eye(numel(tidx));
    else
        tempBlk = Rcorr*ones(cls(it)) + eye(cls(it))*(1-Rcorr);
        [U,S] = eig(tempBlk);
        Blk.U = U;
        Blk.S = diag(S);
        clusterBlk{it} = Blk;
        G(tidx,clusters{it}) = U*diag(sqrt(Blk.S));
    end
end

naSource = nSource+nSensor;
numW = size(Phi,1); % number of TBFs
C_w = cell(numW,1);
C_phi = eye(numW)*1e-6;
alpha = ones(numW,1);
diagCk = zeros(numW,1);
W = ones(sum(cls),numW);
Cgamma = (ones(naSource,1)+1e-3*randn(naSource,1))*trace(LF*LF.')*trace(Phi*Phi.')/(trace(B*B'));
dipIdx = 1:naSource;
prune = [1e-6,1e-1];
cost = -1;

costlist = [];


gamma = [];Ncls = zeros(numel(dipIdx),1);
for it = 1:numel(dipIdx)
    Ncls(it) = numel(clusters{it});
    gamma = [gamma;Cgamma(it)*ones(Ncls(it),1)];
end
remainIdx = 1:sum(cls);

for iter = 1:maxIter
    %% check W and TBF
    tbfIdx = abs(1./alpha)>(max(abs(1./alpha))*prune(2));
    if ~prod(tbfIdx)
        Phi = Phi(tbfIdx,:);
        alpha = alpha(tbfIdx);
        diagCk = diagCk(tbfIdx);
        C_w = C_w(tbfIdx);
        C_phi = C_phi(tbfIdx,tbfIdx);
        numW = sum(tbfIdx);
        W = W(:,tbfIdx);
    end
    
    wIdx = abs(1./Cgamma)>(max(abs(1./Cgamma))*prune(1));
    if ~prod(wIdx)
        dipIdx = dipIdx(wIdx);
        clusters = clusters(wIdx);
        clusterBlk = clusterBlk(wIdx);
        Cgamma = Cgamma(wIdx);
        Ncls = zeros(numel(dipIdx),1);
        gamma = []; remainIdx = [];
        for it = 1:numel(dipIdx)
            Ncls(it) = numel(clusters{it});
            gamma = [gamma;Cgamma(it)*ones(Ncls(it),1)];
            remainIdx = [remainIdx,clusters{it}];
        end
    end
    W = W(remainIdx,:);
    EF = LF*G(:,remainIdx);
    %% spatial coefficient
    FAF = EF.*repmat(1./gamma',nSensor,1)* EF.';
    
    for k = 1:numW
        tIdx = setdiff(1:numW,k);
        x = ((Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k));
        C_w{k} = FAF+C_noise/x;
        
        r_k = B*Phi(k,:).'-EF*sum(W(:,tIdx).*repmat(Phi(k,:)*Phi(tIdx,:).'+nSamp*C_phi(k,tIdx),size(EF,2),1),2);
        
        W(:,k) = repmat(1./gamma,1,nSensor).*EF.'/C_w{k}*r_k/x;
        diagCk(k) = trace(FAF/C_noise-FAF/C_w{k}*FAF/C_noise);% trace(L'L Cwk)
    end
    %% Phi
    LW = EF*W;
    CR = LW.'/C_noise*LW;
    C_phi = inv(CR+diag(alpha+diagCk));
    Phi = C_phi*LW.'/C_noise*(B);
    %% gamma
    mu = zeros(numel(Cgamma),1);
    for it = 1:numel(clusters)
        if it == 1
            col = 1:sum(Ncls(1:it));
        else
            col = sum(Ncls(1:it-1))+1 : sum(Ncls(1:it));
        end
        for k = 1:numW
            x = (Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k);
            C_wk = FAF+C_noise/x; % re-calculate C_w{k}
            mu(it) = mu(it)+trace(EF(:,col)'/C_wk*EF(:,col));
        end
        Cgamma(it) = sqrt(mu(it)/trace(W(col,:)'*W(col,:)));
    end
    %% alpha
    alpha = 1./(diag(C_phi)+diag(Phi*Phi')/nSamp);
    %% free energy
    gamma = [];
    for it = 1:numel(dipIdx)
        gamma = [gamma;Cgamma(it)*ones(Ncls(it),1)];
    end
    tempk = 0;
    temp = 0;
    FAF = EF.*repmat(1./gamma',nSensor,1)*EF.';
    for k = 1:numW
        x = (Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k);
        C_wk = FAF+C_noise/x;
        cws = 1; % prevent Inf results of log(det(C_wk))
        if trace(C_wk)>1e6
            cws = max(10^ceil(log10(trace(C_wk))-6),1e3);
        end
        temp = temp+W(:,k).'.*gamma.'*W(:,k);
        tempk = tempk+log(det(C_wk/cws))+size(C_wk,1)*log(cws)+log(x)*nSensor;
    end
    cost_old = cost;
    cns = 1;
    if trace(C_noise)>1e6
        cns =  max(10^ceil(log10(trace(C_noise))-6),1e3);
    end
    
    cost = -trace((B-LW*Phi)'/C_noise*(B-LW*Phi))-temp-tempk-nSamp*trace(CR*C_phi)...
        +nSamp*(sum(log(alpha))+log(det((C_phi))))-nSamp*(log(det(C_noise*2*pi/cns))+ size(C_noise,1)*log(cns));  
    %% residual noise
    res = sum((B-LW*Phi).^2,2)+diag(LW*C_phi*LW'*nSamp);
    C_noise = diag(max(beta,res/nSamp));    
    costlist = [costlist,cost];     
    %% Recover W
    temp = zeros(sum(cls),numW);
    temp(remainIdx,:) = W;
    W = temp;
    %% the simplified version uses fixed block correlation    
    if abs((cost-cost_old)/cost)<epsilon
        break;
    end
end
%% source
gtW = G*W;
par.W = gtW;
par.Phi = Phi;
par.dipIdx = dipIdx;
par.alpha = alpha;
par.gamma = Cgamma;
par.noise = C_noise;
par.cost = costlist;
S = par.W*Phi;

end


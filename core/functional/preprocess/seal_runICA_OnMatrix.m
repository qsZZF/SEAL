function ica = seal_runICA_OnMatrix( X, alg, rankX, extraArgs)
            alg = validatestring(alg, {'runica','binica','fastica'}, mfilename, 'alg');
            switch alg
                case 'fastica'
                    args = [{'numOfIC', rankX, 'verbose', 'on'}, extraArgs];
                    [S, A, W] = fastica(X, args{:});
                    ica.S = S; ica.W = W; ica.A = A;

                case {'runica','binica'}
                    if strcmp(alg,'binica') && ~exist('binica','file')
                        warning('未找到 binica，改用 runica');
                        alg='runica';
                    end
                    fun = str2func(alg);
                    args = [{'pca', rankX}, extraArgs];
                    argsIn = {};
                    for i = 1:numel(args)
                        v = args{i};
                        while iscell(v) && numel(v)==1    % 递归去壳
                            v = v{1};
                        end
                        argsIn{end+1} = v;
                    end
                    [weights, sphere] = fun(X, argsIn{:});
                    W = weights*sphere;
                    S = W*X;
                    A = pinv(W);
                    ica.S = S; ica.W = W; ica.A = A;
            end
        end
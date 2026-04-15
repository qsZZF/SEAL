function app = SEAL()
%SEAL  Launch SEAL Toolbox GUI
%   在 MATLAB 命令行输入 SEAL 即可启动主界面,模仿 EEGLAB 的启动方式。

    % 1. 把 SEAL 根目录加入 path(保证无论在哪个目录都能找到依赖)
    sealRoot = fileparts(mfilename('fullpath'));
    addpath(sealRoot);
    addpath(genpath(fullfile(sealRoot, 'algorithms')));
    addpath(genpath(fullfile(sealRoot, 'core')));
    addpath(genpath(fullfile(sealRoot, 'external')));
    addpath(genpath(fullfile(sealRoot, 'docImg')));
    addpath(genpath(fullfile(sealRoot, 'docs')));
    addpath(genpath(fullfile(sealRoot, 'utilities')));
    addpath(genpath(fullfile(sealRoot, 'app')));

    % 2. 防重复启动:如果已经开着一个 SEAL 主窗,就把它拉到前台
    existing = findall(0, 'Type', 'figure', 'Name', 'SEAL');
    if ~isempty(existing)
        figure(existing(1));
        if nargout > 0
            app = [];  % 无法拿到已有实例就返回空,或另存一份 handle
        end
        disp('SEAL 已经在运行,已切换到前台。');
        return;
    end

    % 3. 正常启动
    fprintf('Launching SEAL Toolbox...\n');
    appHandle = SEAL_GUI();

    % 4. 只有用户显式接收返回值才输出 app(保持 `SEAL` 裸输入时命令行干净)
    if nargout > 0
        app = appHandle;
    end
end
function chanlocs = ensureFullChanlocs(chanlocs)
%ENSUREFULLCHANLOCS 补全 EEGLAB 标准 chanlocs 字段,避免触发 readlocs
    
    n = numel(chanlocs);
    requiredFields = {'labels', 'theta', 'radius', 'X', 'Y', 'Z', ...
                      'sph_theta', 'sph_phi', 'sph_radius', 'type', 'ref'};
    
    for i = 1:n
        % 补充缺失的字段为空字符串/空数组
        for f = requiredFields
            if ~isfield(chanlocs, f{1}) || isempty(chanlocs(i).(f{1}))
                if strcmp(f{1}, 'labels')
                    chanlocs(i).labels = sprintf('Ch%d', i);
                elseif ismember(f{1}, {'type', 'ref'})
                    chanlocs(i).(f{1}) = '';
                else
                    chanlocs(i).(f{1}) = [];
                end
            end
        end
        
        % 如果只有 X/Y/Z,自动算出 theta/radius
        if isempty(chanlocs(i).theta) && ~isempty(chanlocs(i).X)
            x = chanlocs(i).X;
            y = chanlocs(i).Y;
            z = chanlocs(i).Z;
            % EEGLAB 极坐标定义:theta 是头顶俯视下的方位角(度)
            chanlocs(i).theta = atan2(y, x) * 180 / pi;
            % radius 是从头顶到该点的归一化弧长
            r = sqrt(x^2 + y^2 + z^2);
            if r > 0
                chanlocs(i).radius = (pi/2 - asin(z/r)) / pi;
            else
                chanlocs(i).radius = 0;
            end
        end
        
        % 如果只有 theta/radius,反算 X/Y/Z
        if isempty(chanlocs(i).X) && ~isempty(chanlocs(i).theta)
            th = chanlocs(i).theta * pi / 180;
            r  = chanlocs(i).radius;
            % 投影到单位球
            phi = pi/2 - r * pi;
            chanlocs(i).X = cos(phi) * cos(th);
            chanlocs(i).Y = cos(phi) * sin(th);
            chanlocs(i).Z = sin(phi);
        end
        
        % 球面坐标(可选)
        if isempty(chanlocs(i).sph_theta) && ~isempty(chanlocs(i).X)
            chanlocs(i).sph_theta  = chanlocs(i).theta;
            chanlocs(i).sph_radius = sqrt(chanlocs(i).X^2 + chanlocs(i).Y^2 + chanlocs(i).Z^2);
            chanlocs(i).sph_phi    = asin(chanlocs(i).Z / chanlocs(i).sph_radius) * 180 / pi;
        end
    end
end
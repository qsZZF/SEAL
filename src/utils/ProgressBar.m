classdef ProgressBar < handle
    % MATLAB实现的类似Python tqdm的进度条类
    % 使用动量法计算速度
    
    properties
        total           % 总迭代次数
        current = 0     % 当前迭代次数
        start_time      % 开始时间
        waitbar_handle  % waitbar句柄
        description     % 进度条描述
        unit = 'it'     % 单位
        unit_scale = false % 是否自动缩放单位
        position = []   % 进度条位置
        leave = false   % 完成后是否保留进度条
        auto_close = true % 完成后自动关闭
    end
    
    properties (Access = private)
        last_update_time = 0
        update_interval = 0.1 % 更新间隔(秒)
        
        % 动量法相关参数
        beta = 0.9           % 动量系数，控制历史权重
        smoothed_speed = 0    % 平滑后的速度
        smoothed_speed_initialized = false % 标记是否已初始化平滑速度
        
        last_iteration_time = 0 % 上一次迭代的时间
        is_finished = false   % 标记是否已完成
        user_cancelled = false % 标记是否被用户取消
    end
    
    methods
        function obj = ProgressBar(total, varargin)
            % ProgressBar 构造函数
            obj.total = total;
            obj.start_time = tic;
            obj.last_iteration_time = toc(obj.start_time);
            
            % 解析可选参数
            p = inputParser;
            addParameter(p, 'desc', '', @ischar);
            addParameter(p, 'unit', 'it', @ischar);
            addParameter(p, 'unit_scale', false, @islogical);
            addParameter(p, 'position', [], @(x) isnumeric(x) && length(x) == 2);
            addParameter(p, 'leave', false, @islogical);
            addParameter(p, 'auto_close', true, @islogical);
            addParameter(p, 'beta', 0.9, @(x) isnumeric(x) && x > 0 && x < 1); % 动量系数
            parse(p, varargin{:});
            
            obj.description = p.Results.desc;
            obj.unit = p.Results.unit;
            obj.unit_scale = p.Results.unit_scale;
            obj.position = p.Results.position;
            obj.leave = p.Results.leave;
            obj.auto_close = p.Results.auto_close;
            obj.beta = p.Results.beta;
            
            % 创建初始进度条
            obj.create_waitbar();
        end
        
        function create_waitbar(obj)
            % 创建waitbar
            initial_message = obj.get_message();
            if isempty(obj.position)
                obj.waitbar_handle = waitbar(0, initial_message, ...
                    'Name', 'Progress', 'CreateCancelBtn', @(src,evt) obj.cancel_callback());
            else
                obj.waitbar_handle = waitbar(0, initial_message, ...
                    'Name', 'Progress', 'CreateCancelBtn', @(src,evt) obj.cancel_callback(), ...
                    'Position', [obj.position(1), obj.position(2), 270, 60]);
            end
            
            set(obj.waitbar_handle, 'CloseRequestFcn', @(src,evt) obj.close_request_callback());
        end
        
        function update(obj, n, description)
            % UPDATE 使用动量法计算速度更新进度条
            
            if obj.is_finished || obj.user_cancelled
                return;
            end
            
            if nargin < 2 || isempty(n)
                n = 1;
            end
            
            if nargin < 3
                description = obj.description;
            else
                obj.description = description;
            end
            
            obj.current = obj.current + n;
            
            % 检查是否完成
            if obj.current >= obj.total
                obj.current = obj.total;
                obj.is_finished = true;
            end
            
            % 使用动量法更新速度估计
            current_time = toc(obj.start_time);
            time_diff = current_time - obj.last_iteration_time;
            
            if time_diff > 0
                current_speed = n / time_diff; % 当前瞬时速度
                
                if ~obj.smoothed_speed_initialized
                    % 第一次更新，直接使用当前速度
                    obj.smoothed_speed = current_speed;
                    obj.smoothed_speed_initialized = true;
                else
                    % 使用指数加权移动平均更新速度
                    obj.smoothed_speed = obj.beta * obj.smoothed_speed + (1 - obj.beta) * current_speed;
                end
            end
            
            obj.last_iteration_time = current_time;
            
            % 限制更新频率，这是因为根据实测，如果有快速迭代但是范围很大的循环，反复刷新会严重拖慢性能
            if current_time - obj.last_update_time < obj.update_interval && ~obj.is_finished
                return;
            end
            
            obj.last_update_time = current_time;
            
            % 更新waitbar
            if ishandle(obj.waitbar_handle)
                progress = obj.current / obj.total;
                message = obj.get_message();
                waitbar(progress, obj.waitbar_handle, message);
                
                % 如果完成且设置了自动关闭，则关闭进度条
                if obj.is_finished && obj.auto_close
                    obj.safe_close();
                end
            else
                % 如果句柄无效，标记为已完成
                obj.is_finished = true;
            end
        end
        
        function close(obj)
            % CLOSE 关闭进度条
            obj.safe_close();
        end
        
        function safe_close(obj)
            % SAFE_CLOSE 安全关闭进度条
            if ishandle(obj.waitbar_handle)
                set(obj.waitbar_handle, 'CloseRequestFcn', 'closereq');
                close(obj.waitbar_handle);
            end
            obj.is_finished = true;
            obj.waitbar_handle = [];
        end
        
        function wait_for_close(obj)
            % WAIT_FOR_CLOSE 等待用户手动关闭进度条
            if ishandle(obj.waitbar_handle) && obj.is_finished
                uiwait(obj.waitbar_handle);
            end
        end
        
        function reset(obj, total)
            % RESET 重置进度条
            if nargin > 1
                obj.total = total;
            end
            obj.current = 0;
            obj.start_time = tic;
            obj.last_iteration_time = toc(obj.start_time);
            obj.last_update_time = 0;
            
            % 重置动量法相关状态
            obj.smoothed_speed = 0;
            obj.smoothed_speed_initialized = false;
            
            obj.is_finished = false;
            obj.user_cancelled = false;
            
            if ishandle(obj.waitbar_handle)
                waitbar(0, obj.waitbar_handle, obj.get_message());
            else
                obj.create_waitbar();
            end
        end
        
        function set_description(obj, description)
            % SET_DESCRIPTION 设置描述
            obj.description = description;
            if ishandle(obj.waitbar_handle) && ~obj.is_finished && ~obj.user_cancelled
                waitbar(obj.current/obj.total, obj.waitbar_handle, obj.get_message());
            end
        end
        
        function finish(obj)
            % FINISH 立即完成进度条
            obj.current = obj.total;
            obj.is_finished = true;
            if ishandle(obj.waitbar_handle)
                waitbar(1, obj.waitbar_handle, obj.get_message());
                if obj.auto_close
                    obj.safe_close();
                end
            end
        end
        
        function cancelled = is_cancelled(obj)
            % IS_CANCELLED 检查是否被用户取消
            cancelled = obj.user_cancelled;
        end
        
        function speed = get_smoothed_speed(obj)
            % GET_SMOOTHED_SPEED 获取当前平滑后的速度
            % 这个函数可以用于外部监控进度速度
            if obj.smoothed_speed_initialized
                speed = obj.smoothed_speed;
            else
                speed = 0;
            end
        end
    end
    
    methods (Access = private)
        function msg = get_message(obj)
            % 生成进度条消息
            elapsed = toc(obj.start_time);
            progress = obj.current / obj.total;
            
            % 使用平滑速度计算剩余时间
            if obj.smoothed_speed_initialized && obj.smoothed_speed > 0 && ~obj.is_finished && ~obj.user_cancelled
                remaining = (obj.total - obj.current) / obj.smoothed_speed;
            else
                remaining = 0;
            end
            
            % 格式化输出
            if ~isempty(obj.description)
                desc_str = [obj.description, ': '];
            else
                desc_str = '';
            end
            
            percentage = sprintf('%.1f%%', progress * 100);
            progress_str = sprintf('%d/%d', obj.current, obj.total);
            
            % 格式化时间
            elapsed_str = obj.format_time(elapsed);
            if obj.is_finished
                time_info = sprintf('Total: %s', elapsed_str);
            elseif obj.user_cancelled
                time_info = sprintf('Cancelled: %s', elapsed_str);
            else
                remaining_str = obj.format_time(remaining);
                time_info = sprintf('%s<%s', elapsed_str, remaining_str);
            end
            
            % 格式化速度
            if obj.unit_scale
                rate_str = obj.format_unit(obj.smoothed_speed);
            else
                rate_str = sprintf('%.2f %s/s', obj.smoothed_speed, obj.unit);
            end
            
            % 组合消息
            if obj.user_cancelled
                msg = sprintf('%s Cancelled! %s [%s, %s]', ...
                    desc_str, percentage, progress_str, time_info);
            elseif obj.is_finished
                msg = sprintf('%s Finished! %s [%s, %s]', ...
                    desc_str, percentage, progress_str, time_info);
            else
                msg = sprintf('%s %s %s [%s, %s]', ...
                    desc_str, percentage, progress_str, time_info, rate_str);
            end
        end
        
        function time_str = format_time(~, seconds)
            % 格式化时间显示
            if seconds < 3600
                minutes = floor(seconds / 60);
                seconds = mod(seconds, 60);
                time_str = sprintf('%02d:%02d', minutes, round(seconds));
            else
                hours = floor(seconds / 3600);
                minutes = floor(mod(seconds, 3600) / 60);
                time_str = sprintf('%02d:%02d:%02d', hours, minutes, round(mod(seconds, 60)));
            end
        end
        
        function scaled_str = format_unit(~, value)
            % 自动缩放单位显示
            units = {'', 'K', 'M', 'G', 'T'};
            unit_index = 1;
            while value >= 1000 && unit_index < length(units)
                value = value / 1000;
                unit_index = unit_index + 1;
            end
            scaled_str = sprintf('%.2f %s', value, units{unit_index});
        end
        
        function cancel_callback(obj)
            % 取消按钮回调函数
            obj.user_cancelled = true;
            obj.safe_close();
        end
        
        function close_request_callback(obj)
            % 关闭请求回调函数
            obj.user_cancelled = true;
            obj.safe_close();
        end
    end
    
    methods (Static)
        function t = trange(total, varargin)
            % TRANGE 创建进度条并返回迭代器
            t = ProgressBar(total, varargin{:});
        end
    end
    
    % 析构函数
    methods
        function delete(obj)
            obj.safe_close();
        end
    end
end
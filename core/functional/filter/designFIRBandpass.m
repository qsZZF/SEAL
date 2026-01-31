   function firFilter = designFIRBandpass(Fs, Flow,Fhigh, Rs, N)
            % 设计FIR低通滤波器
            % 输入:
            %   Fs - 采样率 (Hz)
            %   Flow - 低频截止频率 (Hz)
            %   Fhigh - 高频截止频率 (Hz)
            %   Rp - 通带波纹 (dB)
            %   Rs - 阻带衰减 (dB)
            %   N - 滤波器阶数 (必须为偶数)
            % 输出:
            %   firFilter - FIR滤波器系数

            % 参数验证
            assert(Flow > 0, '低频截止频率必须大于0');
            assert(Fhigh > Flow, '高频截止频率必须大于低频截止频率');
            assert(Fhigh < Fs/2, '高频截止频率必须小于奈奎斯特频率(%.1f Hz)', Fs/2);
            assert(mod(N,2) == 0, '阶数N必须是偶数');

            % 归一化频率
            normFlow = Flow / (Fs/2);   % 归一化低频截止频率
            normFhigh = Fhigh / (Fs/2); % 归一化高频截止频率

            % 设计滤波器 (使用Kaiser窗方法)
            beta = kaiserbeta(Rs); % 计算Kaiser窗β参数
            firFilter = fir1(N, [normFlow, normFhigh], 'bandpass', kaiser(N+1, beta), 'noscale');

   end
        
       function beta = kaiserbeta(Rs)
            % 根据阻带衰减计算Kaiser窗的β参数
            if Rs > 50
                beta = 0.1102 * (Rs - 8.7);
            elseif Rs >= 21
                beta = 0.5842 * (Rs - 21)^0.4 + 0.07886 * (Rs - 21);
            else
                beta = 0;
            end
       end
       
       
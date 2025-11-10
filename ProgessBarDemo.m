fprintf('=== ProgressBar使用示例 ===\n');

try
    total = 1000;
    pb1 = ProgressBar(total, 'desc', 'Processing');
    
    for i = 1:total
        pb1.update();

        % 这里只是对计算过程进行模拟，实际替换成需要计算的代码
        % 强烈建议限制刷新频率，此处默认0.1秒刷新一次
        % 理论上每秒可以循环20次，但是在此条件下实际只有17次
        % 这是进度条的开销导致的
        pause(0.05); 
        
        if pb1.is_cancelled()
            fprintf('进度被用户取消\n');
            break;
        end
    end
    
catch ME
    fprintf('发生错误: %s\n', ME.message);
end
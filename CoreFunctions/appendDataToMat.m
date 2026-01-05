function appendDataToMat(protocolpath, fieldName, dataValue)
% 定义保存的文件路径
filename = fullfile(protocolpath, 'ProcessingInfo.mat');

% 1. 确保字段名是合法的变量名 (防止输入空格或特殊字符报错)
validName = matlab.lang.makeValidName(fieldName);


if ~isfile(filename)
    % --- 情况 A：文件不存在 ---
    % 创建新文件并保存
    s = struct();
    s.(validName) = dataValue;
    save(filename, '-struct', 's', '-v7.3');
    return; % 结束函数
end


% 创建 matfile 对象，开启写入权限
m = matfile(filename, 'Writable', true);

% 4. 检查该变量是否已存在于文件中
% who(m) 返回文件内所有变量名的 cell 数组
varList = who(m);

if any(strcmp(varList, validName))
    % === 核心逻辑：变量已存在，进行“扩容” ===

    % 4.1 读取旧数据
    oldData = m.(validName);

    % 4.2 根据类型进行合并
    if isstruct(oldData)
        % --- 处理结构体 (Struct) ---
        % 要求：新旧结构体的字段必须一致，否则会报错
        % 这里的逻辑是将新结构体拼接到旧结构体数组的末尾
        % 自动处理横向或纵向拼接，通常使用 [old, new]
        try
            newData = [oldData, dataValue];
        catch
            % 如果直接拼接失败（通常是因为维度不匹配，比如一个是行向量一个是列向量）
            % 尝试强制转为行向量再拼接
            newData = [oldData(:)', dataValue(:)'];
        end

    elseif iscell(oldData)
        % --- 处理单元数组 (Cell) ---
        % 如果输入的新数据不是 cell，为了放入 cell 数组，通常需要包裹一下
        % 例如：旧的是 {'a'}, 新来的是 'b' -> 变成 {'a', 'b'}
        if ~iscell(dataValue)
            tempValue = {dataValue};
        else
            tempValue = dataValue;
        end

        try
            newData = [oldData, tempValue];
        catch
            % 维度不匹配时，全部转行向量拼接
            newData = [oldData(:)', tempValue(:)'];
        end

    elseif isnumeric(oldData)
        % --- 处理数值矩阵 ---
        % 默认垂直追加 (类似于记录日志)
        newData = [oldData; dataValue];

    else
        % --- 其他类型或混合类型 ---
        % 如果类型不兼容（比如旧的是 double，新的是 char）
        % 强制转为 cell 数组以保存所有数据
        warning('类型不兼容，将转换为 Cell 数组保存。');
        if iscell(oldData)
            newData = [oldData, {dataValue}];
        else
            newData = {oldData, dataValue};
        end
    end

    % 4.3 将合并后的新数据写回文件
    m.(validName) = newData;

else
    % === 变量不存在，直接创建 ===
    m.(validName) = dataValue;
end
end
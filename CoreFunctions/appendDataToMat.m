     function appendDataToMat(protocolpath, fieldName, dataValue)
        % 定义保存的文件路径
        filename = fullfile(protocolpath, 'ProcessingInfo.mat');

        % 创建 matfile 对象，开启写入权限
        m = matfile(filename, 'Writable', true);

        % 1. 确保字段名是合法的变量名 (防止输入空格或特殊字符报错)
        validName = matlab.lang.makeValidName(fieldName);

        % 2. 将数据写入对应的字段
        % 这里的语法 m.(validName) 等同于在 mat 文件中创建了一个名为 validName 的变量
        m.(validName) = dataValue;
        end
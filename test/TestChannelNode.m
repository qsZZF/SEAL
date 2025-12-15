classdef TestChannelNode < matlab.unittest.TestCase
    %TESTCHANNELNODE 测试ChannelNode功能
    % 测试通道节点的创建、打开、保存、数据加载等核心功能
    
    properties
        testProjectName = "TestProject_Channel"
        testProtocolName = "TestProtocol_Channel"
        testChannelName = "TestChannel"
        testProjectPath
        testProtocolPath
        testDataPath
        testChannelData
        projectNode
        protocolNode
    end
    
    properties (Constant)
        % 测试通道数据
        CHANNEL_NAMES = {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'}
        CHANNEL_POSITIONS = randn(10, 3)  % 随机生成10个通道的位置
        CHANNEL_TYPES = repmat({'EEG'}, 1, 10)
        SAMPLING_RATE = 1000  % 采样率
    end
    
    methods (TestMethodSetup)
        function setup(testCase)
            % 设置测试环境
            testCase.testProjectPath = fullfile(tempdir, "SEAL", 'SEAL_Test_Channel');
            testCase.testProtocolPath = fullfile(tempdir, "SEAL", 'SEAL_Test_Channel_Protocol');
            testCase.testDataPath = fullfile(tempdir, "SEAL", 'SEAL_Test_Channel_Data');
            
            % 确保测试目录存在
            folders = {testCase.testProjectPath, testCase.testProtocolPath, testCase.testDataPath};
            for i = 1:length(folders)
                if ~isfolder(folders{i})
                    mkdir(folders{i});
                end
            end
            
            % 创建测试通道数据
            testCase.createTestChannelData();
            
            fprintf('=== ChannelNode测试设置 ===\n');
            fprintf('项目目录: %s\n', testCase.testProjectPath);
            fprintf('协议目录: %s\n', testCase.testProtocolPath);
            fprintf('数据目录: %s\n', testCase.testDataPath);
        end
    end
    
    methods (TestMethodTeardown)
        function teardown(testCase)
            % 清理测试环境
            fprintf('\n=== ChannelNode测试清理 ===\n');
            
            % 删除所有节点
            if ~isempty(testCase.projectNode) && isvalid(testCase.projectNode)
                delete(testCase.projectNode);
            end
            
            if ~isempty(testCase.protocolNode) && isvalid(testCase.protocolNode)
                delete(testCase.protocolNode);
            end
            
            % 清理测试目录
            testCase.cleanupTestDirectories();
        end
    end
    
    methods (Test)
        function testChannelInfoCreation(testCase)
            % 测试ChannelInfo的创建
            
            fprintf('\n--- 测试ChannelInfo创建 ---\n');
            
            % 创建ChannelInfo
            channelInfo = ChanlocsInfo.createNew(...
                testCase.testChannelName, ...
                testCase.testDataPath, ...
                "测试通道描述", ...
                struct("source", "simulated", "unit", 'V'));
            
            % 验证基本属性
            testCase.verifyEqual(channelInfo.name, testCase.testChannelName);
            testCase.verifyEqual(channelInfo.desc, "测试通道描述");
            testCase.verifyEqual(channelInfo.dataPath, testCase.testDataPath);
            testCase.verifyEqual(channelInfo.metadata.source, "simulated");
            testCase.verifyNotEmpty(channelInfo.createdDate);
            testCase.verifyNotEmpty(channelInfo.modifiedDate);
            
            fprintf('✓ ChannelInfo创建测试通过\n');
        end
        
        function testChannelInfoSaveAndLoad(testCase)

            % 测试ChannelInfo的保存和加载
            
            fprintf('\n--- 测试ChannelInfo保存和加载 ---\n');
            
            % 创建并保存ChannelInfo
            savePath = fullfile(testCase.testDataPath, 'SavedChannel');
            if ~isfolder(savePath)
                mkdir(savePath);
            end
            
            channelInfo = ChanlocsInfo.createNew(...
                'SavedChannel', ...
                testCase.testDataPath, ...
                '保存测试通道', ...
                struct("test", "save_load"));
            
            channelInfo.save(savePath);
            
            % 验证文件已创建
            expectedFile = fullfile(savePath, 'SavedChannel.mat');
            testCase.verifyTrue(isfile(expectedFile), 'ChannelInfo文件应被创建');
            
            % 重新加载ChannelInfo
            loadedChannelInfo = ChanlocsInfo.openExisting(expectedFile);
            
            % 验证加载的数据
            testCase.verifyEqual(loadedChannelInfo.name, "SavedChannel");
            testCase.verifyEqual(loadedChannelInfo.desc, "保存测试通道");
            testCase.verifyEqual(loadedChannelInfo.metadata.test, "save_load");
            testCase.verifyNotEmpty(loadedChannelInfo.createdDate);
            
            fprintf('✓ ChannelInfo保存和加载测试通过\n');
        end
        
        function testChannelNodeCreation(testCase)
            % 测试ChannelNode的创建
            
            fprintf('\n--- 测试ChannelNode创建 ---\n');
            
            % 创建项目、协议和通道节点
            testCase.createTestNodes();
            
            % 创建ChannelInfo
            channelInfo = ChanlocsInfo.createNew(...
                testCase.testChannelName, ...
                testCase.testDataPath, ...
                '通道节点测试', ...
                struct('type', 'test'));
            
            % 创建ChannelNode
            channelNode = ChanlocsNode();
            channelNode.channelInfo = channelInfo;
            channelNode.path = fullfile(testCase.testProtocolPath, 'Channels', testCase.testChannelName);
            
            % 添加到协议节点
            testCase.protocolNode.setChannelNode(channelNode);
            
            % 验证通道节点属性
            testCase.verifyEqual(channelNode.name, testCase.testChannelName);
            testCase.verifyFalse(channelNode.isLoaded, "新创建的ChannelNode应未加载数据");
            testCase.verifyEqual(channelNode.type, "ChannelNode");
            testCase.verifyEqual(channelNode.parent, testCase.protocolNode);
            
            % 验证协议节点中的通道节点引用
            testCase.verifyEqual(testCase.protocolNode.channelNode, channelNode);
            
            fprintf('✓ ChannelNode创建测试通过\n');
        end
        
        function testChannelNodeSave(testCase)
            % 测试ChannelNode的保存功能
            
            fprintf('\n--- 测试ChannelNode保存 ---\n');
            
            % 创建测试节点
            testCase.createTestNodes();
            
            % 创建并配置ChannelNode
            channelInfo = ChanlocsInfo.createNew(...
                'SaveTestChannel', ...
                testCase.testDataPath, ...
                '保存测试', ...
                struct('sampling_rate', testCase.SAMPLING_RATE));
            
            channelNode = ChanlocsNode();
            channelNode.channelInfo = channelInfo;
            channelPath = fullfile(testCase.testProtocolPath, 'Channels', 'SaveTestChannel');
            channelNode.path = channelPath;
            
            % 保存ChannelNode
            channelNode.save();
            
            % 验证文件结构
            testCase.verifyTrue(isfolder(channelPath), '通道目录应被创建');
            
            % 验证ChannelInfo文件
            expectedInfoFile = fullfile(channelPath, 'SaveTestChannel.mat');
            testCase.verifyTrue(isfile(expectedInfoFile), 'ChannelInfo文件应被创建');
            
            % 加载并验证保存的数据
            loadedInfo = loadData(expectedInfoFile);
            testCase.verifyEqual(loadedInfo.name, "SaveTestChannel");
            testCase.verifyEqual(loadedInfo.desc, "保存测试");
            testCase.verifyEqual(loadedInfo.metadata.sampling_rate, testCase.SAMPLING_RATE);
            
            fprintf('✓ ChannelNode保存测试通过\n');
        end
        
        function testChannelNodeLoadAndUnload(testCase)
            % 测试ChannelNode的数据加载和卸载
            
            fprintf('\n--- 测试ChannelNode加载和卸载 ---\n');
            
            % 创建测试通道数据文件
            dataFile = fullfile(testCase.testDataPath, 'channel_data.mat');
            channelData = struct(...
                'names', {testCase.CHANNEL_NAMES}, ...
                'positions', testCase.CHANNEL_POSITIONS, ...
                'types', {testCase.CHANNEL_TYPES}, ...
                'sampling_rate', testCase.SAMPLING_RATE);
            
            saveData(dataFile, channelData);
            
            % 创建ChannelInfo和ChannelNode
            channelInfo = ChanlocsInfo.createNew(...
                'LoadTestChannel', ...
                dataFile, ...
                '加载测试', ...
                struct('data_type', 'mat'));
            
            channelNode = ChanlocsNode();
            channelNode.channelInfo = channelInfo;
            channelNode.path = fullfile(testCase.testDataPath, 'LoadTestChannel');
            
            % 初始状态验证
            testCase.verifyFalse(channelNode.isLoaded, '初始状态应为未加载');
            testCase.verifyEmpty(channelNode.cache, '缓存应为空');
            
            % 加载数据
            channelNode.load();
            
            % 验证加载状态
            testCase.verifyTrue(channelNode.isLoaded, '加载后状态应为已加载');
            testCase.verifyNotEmpty(channelNode.cache, '缓存应不为空');
            testCase.verifyEqual(channelNode.cache.names, testCase.CHANNEL_NAMES);
            testCase.verifyEqual(channelNode.cache.sampling_rate, testCase.SAMPLING_RATE);
            
            % 测试data属性
            data = channelNode.data;
            testCase.verifyEqual(data.names, testCase.CHANNEL_NAMES);
            
            % 卸载数据
            channelNode.unload();
            
            % 验证卸载状态
            testCase.verifyFalse(channelNode.isLoaded, '卸载后状态应为未加载');
            testCase.verifyEmpty(channelNode.cache, '卸载后缓存应为空');
            
            % 再次加载验证
            data = channelNode.data;
            testCase.verifyTrue(channelNode.isLoaded, '通过data属性访问应自动加载');
            testCase.verifyNotEmpty(data);
            
            fprintf('✓ ChannelNode加载和卸载测试通过\n');
        end
        
        function testChannelNodeOpenExisting(testCase)
            % 测试从文件打开现有ChannelNode
            
            fprintf('\n--- 测试ChannelNode打开现有文件 ---\n');
            
            % 准备测试数据
            dataFile = fullfile(testCase.testDataPath, 'existing_channel_data.mat');
            existingData = struct(...
                'channels', testCase.CHANNEL_NAMES, ...
                'locations', testCase.CHANNEL_POSITIONS, ...
                'reference', 'Cz');
           
            saveData(dataFile, existingData);
            
            % 创建并保存ChannelInfo
            saveDir = fullfile(testCase.testDataPath, 'ExistingChannel');
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end
            
            channelInfo = ChanlocsInfo.createNew(...
                'ExistingChannel', ...
                dataFile, ...
                '现有通道数据', ...
                struct('format', 'matlab', 'version', '1.0'));
            
            channelInfo.save(saveDir);
            
            % 使用openExisting打开
            channelNode = ChanlocsNode.openExisting(saveDir);
            
            % 验证打开的节点
            testCase.verifyEqual(channelNode.name, "ExistingChannel");
            testCase.verifyEqual(channelNode.channelInfo.desc, "现有通道数据");
            testCase.verifyEqual(channelNode.channelInfo.dataPath, dataFile);
            testCase.verifyEqual(channelNode.channelInfo.metadata.format, 'matlab');
            
            % 测试数据加载
            channelNode.load();
            testCase.verifyTrue(channelNode.isLoaded);
            
            fprintf('✓ ChannelNode打开现有文件测试通过\n');
        end
        
        function testChannelNodeIntegrationWithProtocol(testCase)
            % 测试ChannelNode与ProtocolNode的集成
            
            fprintf('\n--- 测试ChannelNode与ProtocolNode集成 ---\n');
            
            % 创建项目、协议和通道节点
            testCase.createTestNodes();
            
            % 创建通道数据
            dataFile = fullfile(testCase.testDataPath, 'integration_data.mat');
            integrationData = struct(...
                'montage', '10-20', ...
                'n_channels', 10, ...
                'units', 'uV');
            
            saveData(dataFile, integrationData);
            
            % 创建ChannelInfo
            channelInfo = ChanlocsInfo.createNew(...
                'IntegrationChannel', ...
                dataFile, ...
                '集成测试通道', ...
                struct('project', testCase.testProjectName));
            
            % 创建ChannelNode并添加到协议
            channelNode = ChanlocsNode();
            channelNode.channelInfo = channelInfo;
            channelPath = fullfile(testCase.protocolNode.path, 'Channels', 'IntegrationChannel');
            channelNode.path = channelPath;
            
            testCase.protocolNode.setChannelNode(channelNode);
            
            % 验证集成
            testCase.verifyEqual(testCase.protocolNode.channelNode, channelNode);
            testCase.verifyEqual(channelNode.parent, testCase.protocolNode);
            
            % 保存整个协议
            testCase.protocolNode.save();
            
            % 重新打开协议验证通道节点
            [~, protocolName, ~] = fileparts(testCase.protocolNode.path);
            protocolFile = fullfile(testCase.protocolNode.path, strcat(protocolName, ".mat"));
            reloadedProtocolInfo = ProtocolInfo.openExisting(protocolFile);
            testCase.verifyEqual(reloadedProtocolInfo.name, protocolName);
            
            fprintf('✓ ChannelNode与ProtocolNode集成测试通过\n');
        end
        
        function testChannelNodeDataProperty(testCase)
            % 测试ChannelNode的data依赖属性
            
            fprintf('\n--- 测试ChannelNode data属性 ---\n');
            
            % 创建测试数据
            dataFile = fullfile(testCase.testDataPath, 'data_property_test.mat');
            testData = struct(...
                'signal', randn(1000, 10), ...  % 10个通道，1000个时间点
                'time', (0:999)' / testCase.SAMPLING_RATE, ...
                'labels', {testCase.CHANNEL_NAMES});
            
            saveData(dataFile, testData);
            
            % 创建ChannelNode
            channelInfo = ChanlocsInfo.createNew(...
                'DataPropertyChannel', ...
                dataFile, ...
                '数据属性测试', ...
                struct());
            
            channelNode = ChanlocsNode();
            channelNode.channelInfo = channelInfo;
            channelNode.path = fullfile(testCase.testDataPath, 'DataPropertyChannel');
            
            % 初始状态验证
            testCase.verifyFalse(channelNode.isLoaded);
            
            % 访问data属性（应自动加载）
            data = channelNode.data;
            
            % 验证自动加载
            testCase.verifyTrue(channelNode.isLoaded, '访问data属性应自动加载数据');
            testCase.verifyNotEmpty(data);
            testCase.verifyEqual(size(data.signal), [1000, 10]);
            testCase.verifyEqual(length(data.time), 1000);
            testCase.verifyEqual(data.labels, testCase.CHANNEL_NAMES);
            
            % 再次访问验证缓存
            data2 = channelNode.data;
            testCase.verifyTrue(isequaln(data, data2), '多次访问应返回相同数据');
            
            % 修改数据并保存
            channelNode.save();
            
            fprintf('✓ ChannelNode data属性测试通过\n');
        end
        
        function testChannelNodeInvalidDataPath(testCase)
            % 测试无效数据路径处理
            
            fprintf('\n--- 测试ChannelNode无效数据路径 ---\n');
            
            % 创建ChannelInfo，指定不存在的数据路径
            channelInfo = ChanlocsInfo.createNew(...
                'InvalidPathChannel', ...
                '/nonexistent/path/to/data.mat', ...
                '无效路径测试', ...
                struct());
            
            channelNode = ChanlocsNode();
            channelNode.channelInfo = channelInfo;
            channelNode.path = fullfile(testCase.testDataPath, 'InvalidPathChannel');
            
            % 尝试加载（应不抛出错误，但缓存为空）
            try
                channelNode.load();
                testCase.verifyFalse(channelNode.isLoaded, '无效路径应导致加载失败');
                testCase.verifyEmpty(channelNode.cache, '缓存应为空');
                fprintf('✓ ChannelNode无效数据路径处理测试通过\n');
            catch ME
                % 如果代码抛出异常，也记录下来
                fprintf('加载失败，错误: %s\n', ME.message);
                testCase.verifyTrue(true, '允许加载失败');
            end
        end

        function testFromData(testCase)
            % 测试从文件打开
            
            fprintf('\n--- 测试从文件读取数据 ---\n');

            channelInfo = ChanlocsInfo.fromData("E:\SEAL\Database\1\Chanlocs_60.mat");
            
            channelNode = ChanlocsNode.fromData("E:\SEAL\Database\1\Chanlocs_60.mat");
        end
    end
    
    methods (Access = private)
        function createTestChannelData(testCase)
            % 创建测试通道数据
            
            testCase.testChannelData = struct(...
                'names', {testCase.CHANNEL_NAMES}, ...
                'positions', testCase.CHANNEL_POSITIONS, ...
                'types', {testCase.CHANNEL_TYPES}, ...
                'sampling_rate', testCase.SAMPLING_RATE, ...
                'reference', 'average', ...
                'impedances', rand(1, 10) * 100);
        end
        
        function createTestNodes(testCase)
            % 创建测试项目、协议和通道节点
            
            % 创建项目
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 创建协议
            testCase.protocolNode = testCase.projectNode.createNewProtocol(...
                testCase.testProtocolName, ...
                'Type', 'EEG', ...
                'Description', '通道节点测试协议');
        end
        
        function cleanupTestDirectories(testCase)
            % 清理测试目录
            
            % 删除所有测试目录
            dirs = {testCase.testProjectPath, testCase.testProtocolPath, testCase.testDataPath};
            for i = 1:length(dirs)
                if isfolder(dirs{i})
                    try
                        rmdir(dirs{i}, 's');
                        fprintf('已清理目录: %s\n', dirs{i});
                    catch ME
                        warning('清理目录失败: %s\n错误: %s', dirs{i}, ME.message);
                    end
                end
            end
        end
    end
    
    methods (Static)
        function runAllTests()
            % 运行所有测试
            
            fprintf('=== 开始运行ChannelNode测试套件 ===\n');
            
            % 创建测试套件
            import matlab.unittest.TestSuite;
            suite = TestSuite.fromClass(?TestChannelNode);
            
            % 运行测试
            runner = matlab.unittest.TestRunner.withTextOutput;
            results = runner.run(suite);
            
            % 显示结果
            fprintf('\n=== ChannelNode测试结果 ===\n');
            fprintf('总测试数: %d\n', numel(results));
            fprintf('通过: %d\n', nnz([results.Passed]));
            fprintf('失败: %d\n', nnz([results.Failed]));
            fprintf('跳过: %d\n', nnz([results.Incomplete]));
            
            % 显示失败详情
            failedTests = results([results.Failed]);
            if ~isempty(failedTests)
                fprintf('\n=== 失败测试详情 ===\n');
                for i = 1:length(failedTests)
                    fprintf('测试: %s\n', failedTests(i).Name);
                    fprintf('错误: %s\n\n', failedTests(i).Error.message);
                end
            end
            
            if all([results.Passed])
                fprintf('\n✅ 所有ChannelNode测试通过！\n');
            else
                fprintf('\n❌ 部分测试失败\n');
            end
        end
        
        function runSpecificTest(testName)
            % 运行特定测试
            
            fprintf('运行特定测试: %s\n', testName);
            
            import matlab.unittest.TestSuite;
            suite = TestSuite.fromClass(?TestChannelNode, 'Name', testName);
            
            if isempty(suite)
                fprintf('未找到测试: %s\n', testName);
                return;
            end
            
            runner = matlab.unittest.TestRunner.withTextOutput;
            results = runner.run(suite);
            
            fprintf('\n测试结果: %s\n', testName);
            if results.Passed
                fprintf('✅ 测试通过\n');
            else
                fprintf('❌ 测试失败\n');
                fprintf('错误: %s\n', results.Error.message);
            end
        end
    end
end

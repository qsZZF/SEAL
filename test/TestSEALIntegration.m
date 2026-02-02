classdef TestSEALIntegration < matlab.unittest.TestCase
    %TESTSEALINTEGRATION SEAL框架集成测试
    % 测试整个SEAL框架的数据节点（Project、Protocol、Session、Chanlocs、Cortex、Leadfield、Data）
    % 的创建、读取、保存等集成功能
    
    properties
        testRootPath
        testProjectName = "TestIntegrationProject"
        testProtocolName = "TestIntegrationProtocol"
        testSessionName = "TestIntegrationSession"
        testDataNames = ["original_data", "processed_data_1", "processed_data_2"]
        
        testChanlocsName = "TestChanlocss"
        testCortexName = "TestCortex"
        testLeadfieldName = "TestLeadfield"
        
        projectNode
        protocolNode
        sessionNode
    end
    
    methods (TestClassSetup)
        function setupClass(testCase)
            % 类级别的设置，创建测试数据目录和模拟数据
            
            % 创建测试根目录
            testCase.testRootPath = fullfile(tempdir, "SEAL", "IntegrationTests");
            
            % 创建模拟数据目录
            testDataPath = fullfile(testCase.testRootPath, "TestData");
            if ~isfolder(testDataPath)
                mkdir(testDataPath);
            end
            
            fprintf("=== SEAL集成测试类初始化 ===\n");
            fprintf("测试根目录: %s\n", testCase.testRootPath);
            fprintf("测试数据目录: %s\n", testDataPath);
            
            % 创建测试通道数据
            chanlocsFilePath = fullfile(testDataPath, "test_chanlocss.mat");
            if ~isfile(chanlocsFilePath)
                fprintf("创建测试通道数据: %s\n", chanlocsFilePath);
                chanlocsData = struct(...
                    "labels", {"Fz", "Cz", "Pz", "Oz", "F3", "F4", "C3", "C4", "P3", "P4"}, ...
                    "positions", rand(10, 3), ...
                    "types", repmat({"EEG"}, 1, 10), ...
                    "metadata", struct("samplingRate", 1000, "unit", "mm", "reference", "Cz"));
                save(chanlocsFilePath, "chanlocsData");
            end
            
            % 创建测试皮层数据
            cortexFilePath = fullfile(testDataPath, "test_cortex.mat");
            if ~isfile(cortexFilePath)
                fprintf("创建测试皮层数据: %s\n", cortexFilePath);
                nVertices = 200;
                nFaces = 396;
                cortexData = struct(...
                    "vertices", rand(nVertices, 3), ...
                    "faces", randi([1, nVertices], nFaces, 3), ...
                    "normals", rand(nVertices, 3), ...
                    "metadata", struct("source", "template", "unit", "mm", "subject", "standard"));
                save(cortexFilePath, "cortexData");
            end
            
            % 创建测试导联场数据
            leadfieldFilePath = fullfile(testDataPath, "test_leadfield.mat");
            if ~isfile(leadfieldFilePath)
                fprintf("创建测试导联场数据: %s\n", leadfieldFilePath);
                leadfieldData = struct(...
                    "data", rand(10, nVertices), ...
                    "chanlocs_names", {"Fz", "Cz", "Pz", "Oz", "F3", "F4", "C3", "C4", "P3", "P4"}, ...
                    "metadata", struct("method", "BEM", "conductivity", [0.33, 0.0042, 0.33], ...
                    "resolution", "1mm", "source_space", "cortex"));
                save(leadfieldFilePath, "leadfieldData");
            end
            
            % 创建测试会话数据（多个）
            for i = 1:length(testCase.testDataNames)
                dataFilePath = fullfile(testDataPath, testCase.testDataNames(i) + ".mat");
                if ~isfile(dataFilePath)
                    fprintf("创建测试数据文件: %s\n", dataFilePath);
                    if testCase.testDataNames(i) == "original_data"
                        % 原始数据 - 模拟EEG数据
                        data = struct(...
                            "eeg_data", randn(10, 1000), ...
                            "time", linspace(0, 1, 1000), ...
                            "metadata", struct("task", "resting", "condition", "eyes_open", ...
                            "duration", 1.0, "sampling_rate", 1000));
                    else
                        % 处理后的数据
                        data = struct(...
                            "processed_data", randn(10, 500), ...
                            "time", linspace(0, 0.5, 500), ...
                            "metadata", struct("processing_step", testCase.testDataNames(i), ...
                            "algorithm", "ICA", "components", 10));
                    end
                    save(dataFilePath, "data");
                end
            end
            
            fprintf("=== 测试数据创建完成 ===\n\n");
        end
    end
    
    methods (TestMethodSetup)
        function setup(testCase)
            % 每个测试方法前的设置
            fprintf("\n=== 开始新测试方法 ===\n");
            
            % 清理之前的测试目录
            testProjectPath = fullfile(testCase.testRootPath, testCase.testProjectName);
            if isfolder(testProjectPath)
                rmdir(testProjectPath, "s");
            end
            
            % 重置节点引用
            testCase.projectNode = [];
            testCase.protocolNode = [];
            testCase.sessionNode = [];
            
            fprintf("测试环境已重置\n");
        end
    end
    
    methods (TestMethodTeardown)
        function teardown(testCase)
            % 每个测试方法后的清理
            if ~isempty(testCase.projectNode) && isvalid(testCase.projectNode)
                testCase.projectNode.unload();
                delete(testCase.projectNode);
            end
            
            fprintf("测试完成，清理环境\n");
        end
    end
    
    methods (TestClassTeardown)
        function teardownClass(testCase)
            % 类级别的清理
            if isfolder(testCase.testRootPath)
                try
                    rmdir(testCase.testRootPath, "s");
                    fprintf("清理测试根目录: %s\n", testCase.testRootPath);
                catch ME
                    warning("清理目录失败: %s");
                end
            end
        end
    end
    
    methods (Test)
        function testCreateFullProjectHierarchy(testCase)
            % 测试创建完整的项目层次结构
            
            fprintf("\n=== 测试创建完整项目层次结构 ===\n");
            
            % 1. 创建项目
            fprintf("1. 创建项目: %s\n", testCase.testProjectName);
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testRootPath, ...
                "Description", "集成测试项目");
            
            % 验证项目创建
            testCase.verifyEqual(testCase.projectNode.name, testCase.testProjectName);
            testCase.verifyTrue(isfolder(testCase.projectNode.path), "项目目录应该存在");
            projectInfoFile = fullfile(testCase.projectNode.path, testCase.testProjectName + ".mat");
            testCase.verifyTrue(isfile(projectInfoFile), "项目信息文件应该存在");
            
            % 2. 创建协议
            fprintf("2. 创建协议: %s\n", testCase.testProtocolName);
            testCase.protocolNode = testCase.projectNode.createNewProtocol(...
                testCase.testProtocolName, ...
                "Type", "EEG_STUDY", ...
                "Description", "集成测试协议");
            
            % 验证协议创建
            testCase.verifyEqual(testCase.protocolNode.name, testCase.testProtocolName);
            testCase.verifyTrue(isfolder(testCase.protocolNode.path), "协议目录应该存在");
            protocolInfoFile = fullfile(testCase.protocolNode.path, testCase.testProtocolName + ".mat");
            testCase.verifyTrue(isfile(protocolInfoFile), "协议信息文件应该存在");
            
            % 3. 为协议添加单例节点（Chanlocs, Cortex, Leadfield）
            testDataPath = fullfile(testCase.testRootPath, "TestData");
            
            fprintf("3.1 添加通道节点\n");
            chanlocsDataPath = fullfile(testDataPath, "test_chanlocss.mat");
            testCase.protocolNode.openChanlocsFromData(chanlocsDataPath);
            testCase.verifyNotEmpty(testCase.protocolNode.chanlocsNode, "通道节点应该存在");
            testCase.verifyEqual(testCase.protocolNode.chanlocsNode.name, "test_chanlocss");
            testCase.verifyTrue(isfile(fullfile(testCase.protocolNode.path, "chanlocs_info.mat")), ...
                "通道信息文件应该存在");
            
            fprintf("3.2 添加皮层节点\n");
            cortexDataPath = fullfile(testDataPath, "test_cortex.mat");
            testCase.protocolNode.openCortexFromData(cortexDataPath);
            testCase.verifyNotEmpty(testCase.protocolNode.cortexNode, "皮层节点应该存在");
            testCase.verifyEqual(testCase.protocolNode.cortexNode.name, "test_cortex");
            testCase.verifyTrue(isfile(fullfile(testCase.protocolNode.path, "cortex_info.mat")), ...
                "皮层信息文件应该存在");
            
            fprintf("3.3 添加导联场节点\n");
            leadfieldDataPath = fullfile(testDataPath, "test_leadfield.mat");
            testCase.protocolNode.openLeadfieldFromData(leadfieldDataPath);
            testCase.verifyNotEmpty(testCase.protocolNode.leadfieldNode, "导联场节点应该存在");
            testCase.verifyEqual(testCase.protocolNode.leadfieldNode.name, "test_leadfield");
            testCase.verifyTrue(isfile(fullfile(testCase.protocolNode.path, "leadfield_info.mat")), ...
                "导联场信息文件应该存在");
            
            % 4. 创建会话
            fprintf("4. 创建会话: %s\n", testCase.testSessionName);
            testCase.sessionNode = SessionNode.createNew(...
                testCase.testSessionName, ...
                fullfile(testCase.protocolNode.path, "Sessions"), ...
                "Type", "EEG_RECORDING", ...
                "Description", "集成测试会话");
            
            % 验证会话创建
            testCase.verifyEqual(testCase.sessionNode.name, testCase.testSessionName);
            testCase.verifyTrue(isfolder(testCase.sessionNode.path), "会话目录应该存在");
            sessionInfoFile = fullfile(testCase.sessionNode.path, "session_info.mat");
            testCase.verifyTrue(isfile(sessionInfoFile), "会话信息文件应该存在");
            
            % 5. 为会话添加数据节点
            fprintf("5. 为会话添加数据节点\n");
            for i = 1:length(testCase.testDataNames)
                dataName = testCase.testDataNames(i);
                dataPath = fullfile(testDataPath, dataName + ".mat");
                
                fprintf("  5.%d 添加数据节点: %s\n", i, dataName);
                testCase.sessionNode.openDataFromData(dataPath);
                
                % 验证数据节点添加
                testCase.verifyEqual(testCase.sessionNode.children(i).name, dataName);
            end
            
            % 验证会话数据节点数量
            testCase.verifyEqual(testCase.sessionNode.childCount, length(testCase.testDataNames));
            
            % 6. 保存整个项目
            fprintf("6. 保存整个项目\n");
            testCase.projectNode.save();
            
            % 验证所有文件都存在
            testCase.verifyProjectFileStructure();
            
            fprintf("✓ 完整项目层次结构创建测试通过\n");
        end
        
        function testLoadAndVerifyFullProject(testCase)
            % 测试加载完整项目并验证数据一致性
            
            fprintf("\n=== 测试加载完整项目并验证数据一致性 ===\n");
            
            % 首先创建完整项目
            fprintf("1. 首先创建完整项目\n");
            testCase.testCreateFullProjectHierarchy();
            
            % 保存项目引用
            originalProject = testCase.projectNode;
            
            % 重置节点引用
            testCase.projectNode = [];
            testCase.protocolNode = [];
            testCase.sessionNode = [];
            
            % 2. 重新加载项目
            fprintf("2. 重新加载项目\n");
            projectPath = originalProject.path;
            testCase.projectNode = ProjectNode.openExisting(projectPath);
            
            % 验证项目加载
            testCase.verifyEqual(testCase.projectNode.name, testCase.testProjectName);
            testCase.verifyEqual(testCase.projectNode.protocolCount, 1);
            
            % 获取协议节点
            testCase.protocolNode = testCase.projectNode.getChild(testCase.testProtocolName);
            testCase.verifyNotEmpty(testCase.protocolNode);
            testCase.verifyEqual(testCase.protocolNode.name, testCase.testProtocolName);
            
            % 验证单例节点加载
            testCase.verifyNotEmpty(testCase.protocolNode.chanlocsNode);
            testCase.verifyNotEmpty(testCase.protocolNode.cortexNode);
            testCase.verifyNotEmpty(testCase.protocolNode.leadfieldNode);
            
            % 验证会话数量
            testCase.verifyEqual(testCase.protocolNode.sessionCount, 1);
            
            % 获取会话节点
            testCase.sessionNode = testCase.protocolNode.getChild(testCase.testSessionName);
            testCase.verifyNotEmpty(testCase.sessionNode);
            testCase.verifyEqual(testCase.sessionNode.name, testCase.testSessionName);
            
            % 验证会话数据节点数量
            testCase.verifyEqual(testCase.sessionNode.childCount, length(testCase.testDataNames));
            
            % 3. 加载数据并验证一致性
            fprintf("3. 加载数据并验证一致性\n");
            
            % 加载协议数据
            testCase.protocolNode.load();
            
            % 验证通道数据加载
            if ~isempty(testCase.protocolNode.chanlocsNode)
                testCase.protocolNode.chanlocsNode.load();
                testCase.verifyTrue(testCase.protocolNode.chanlocsNode.isLoaded);
                testCase.verifyNotEmpty(testCase.protocolNode.chanlocsNode.data);
                
                % 检查通道数据是否包含预期字段
                chanlocsData = testCase.protocolNode.chanlocsNode.data;
                testCase.verifyTrue(isfield(chanlocsData, "labels"));
                testCase.verifyTrue(isfield(chanlocsData, "positions"));
                testCase.verifyTrue(isfield(chanlocsData, "metadata"));
            end
            
            % 验证皮层数据加载
            if ~isempty(testCase.protocolNode.cortexNode)
                testCase.protocolNode.cortexNode.load();
                testCase.verifyTrue(testCase.protocolNode.cortexNode.isLoaded);
                testCase.verifyNotEmpty(testCase.protocolNode.cortexNode.data);
                
                % 检查皮层数据是否包含预期字段
                cortexData = testCase.protocolNode.cortexNode.data;
                testCase.verifyTrue(isfield(cortexData, "vertices"));
                testCase.verifyTrue(isfield(cortexData, "faces"));
                testCase.verifyTrue(isfield(cortexData, "metadata"));
            end
            
            % 验证导联场数据加载
            if ~isempty(testCase.protocolNode.leadfieldNode)
                testCase.protocolNode.leadfieldNode.load();
                testCase.verifyTrue(testCase.protocolNode.leadfieldNode.isLoaded);
                testCase.verifyNotEmpty(testCase.protocolNode.leadfieldNode.data);
                
                % 检查导联场数据是否包含预期字段
                leadfieldData = testCase.protocolNode.leadfieldNode.data;
                testCase.verifyTrue(isfield(leadfieldData, "data"));
                testCase.verifyTrue(isfield(leadfieldData, "metadata"));
            end
            
            % 验证会话数据节点加载
            testCase.sessionNode.load();
            for i = 1:testCase.sessionNode.childCount
                dataNode = testCase.sessionNode.children(i);
                dataNode.load();
                testCase.verifyTrue(dataNode.isLoaded);
                testCase.verifyNotEmpty(dataNode.data);
            end
            
            % 4. 验证数据一致性（比较原始数据和加载的数据）
            fprintf("4. 验证数据一致性\n");
            
            % 重新加载原始数据文件进行比较
            testDataPath = fullfile(testCase.testRootPath, "TestData");
            
            % 验证通道数据一致性
            if ~isempty(testCase.protocolNode.chanlocsNode)
                originalChanlocsData = load(fullfile(testDataPath, "test_chanlocss.mat"));
                loadedChanlocsData = testCase.protocolNode.chanlocsNode.data;
                
                testCase.verifyEqual(length({loadedChanlocsData.labels}), length({originalChanlocsData.chanlocsData.labels}));
            end
            
            fprintf("✓ 加载完整项目并验证数据一致性测试通过\n");
        end
        
        function testAddRemoveNodes(testCase)
            % 测试添加和移除节点
            
            fprintf("\n=== 测试添加和移除节点 ===\n");
            
            % 创建基本项目结构
            testCase.createBasicProjectStructure();
            
            % 1. 测试添加新会话
            fprintf("1. 测试添加新会话\n");
            newSessionName = "AdditionalSession";
            newSessionNode = SessionNode.createNew(...
                newSessionName, ...
                fullfile(testCase.protocolNode.path, "Sessions"), ...
                "Type", "EEG_RECORDING", ...
                "Description", "附加测试会话");
            
            testCase.protocolNode.addChild(newSessionNode);
            testCase.verifyEqual(testCase.protocolNode.sessionCount, 2);
            
            % 2. 测试添加新数据节点到会话
            fprintf("2. 测试添加新数据节点到会话\n");
            testDataPath = fullfile(testCase.testRootPath, "TestData");
            newDataPath = fullfile(testDataPath, "test_chanlocss.mat"); % 重用现有文件
            
            newDataInfo = DataInfo.fromData(newDataPath, "Name", "additional_data");
            newDataNode = DataNode.fromInfo(newDataInfo);
            newDataNode.path = newSessionNode.path;
            newSessionNode.addChild(newDataNode);
            
            testCase.verifyEqual(newSessionNode.childCount, 1);
            testCase.verifyEqual(newSessionNode.children(1).name, "additional_data");
            
            % 3. 测试移除数据节点
            fprintf("3. 测试移除数据节点\n");
            originalCount = newSessionNode.childCount;
            newSessionNode.removeChild(newDataNode);
            testCase.verifyEqual(newSessionNode.childCount, originalCount - 1);
            
            % 4. 测试移除会话
            fprintf("4. 测试移除会话\n");
            originalProtocolChildCount = testCase.protocolNode.childCount;
            testCase.protocolNode.removeChild(newSessionNode);
            testCase.verifyEqual(testCase.protocolNode.childCount, originalProtocolChildCount - 1);
            
            % 5. 测试添加第二个协议
            fprintf("5. 测试添加第二个协议\n");
            secondProtocolName = "SecondProtocol";
            secondProtocolNode = testCase.projectNode.createNewProtocol(...
                secondProtocolName, ...
                "Type", "MEG_STUDY", ...
                "Description", "第二个测试协议");
            
            testCase.verifyEqual(testCase.projectNode.protocolCount, 2);
            testCase.verifyEqual(testCase.projectNode.children(2).name, secondProtocolName);
            
            % 保存项目
            testCase.projectNode.save();
            
            fprintf("✓ 添加和移除节点测试通过\n");
        end
        
        function testSaveAndLoadWithModifications(testCase)
            % 测试修改后保存并重新加载
            
            fprintf("\n=== 测试修改后保存并重新加载 ===\n");
            
            % 创建基本项目结构
            testCase.createBasicProjectStructure();
            
            % 1. 修改项目信息
            fprintf("1. 修改项目信息\n");
            testCase.projectNode.projectInfo.desc = "修改后的项目描述";
            
            % 2. 修改协议信息
            fprintf("2. 修改协议信息\n");
            testCase.protocolNode.protocolInfo.desc = "修改后的协议描述";
            testCase.protocolNode.protocolInfo.type = "MEG_EEG_STUDY";
            
            % 3. 修改会话信息
            fprintf("3. 修改会话信息\n");
            testCase.sessionNode.sessionInfo.desc = "修改后的会话描述";
            testCase.sessionNode.sessionInfo.type = "MEG_RECORDING";
            
            % 4. 添加一些自定义元数据
            fprintf("4. 添加自定义元数据\n");
            testCase.projectNode.projectInfo.metadata.testField = "testValue";
            testCase.projectNode.projectInfo.metadata.timestamp = datetime("now");
            
            testCase.protocolNode.protocolInfo.metadata.experimenter = "TestUser";
            testCase.protocolNode.protocolInfo.metadata.version = "1.0.0";
            
            % 5. 保存修改
            fprintf("5. 保存修改\n");
            testCase.projectNode.save();
            
            % 保存项目路径和名称用于重新加载
            projectPath = testCase.projectNode.path;
            
            % 重置节点引用
            testCase.projectNode = [];
            testCase.protocolNode = [];
            testCase.sessionNode = [];
            
            % 6. 重新加载项目
            fprintf("6. 重新加载项目\n");
            testCase.projectNode = ProjectNode.openExisting(projectPath);
            
            % 获取协议节点
            testCase.protocolNode = testCase.projectNode.getChild(testCase.testProtocolName);
            testCase.verifyNotEmpty(testCase.protocolNode);
            
            % 验证协议修改
            testCase.verifyEqual(testCase.protocolNode.protocolInfo.desc, "修改后的协议描述");
            testCase.verifyEqual(testCase.protocolNode.protocolInfo.type, "MEG_EEG_STUDY");
            
            % 验证协议元数据
            testCase.verifyTrue(isfield(testCase.protocolNode.protocolInfo.metadata, "experimenter"));
            testCase.verifyEqual(testCase.protocolNode.protocolInfo.metadata.experimenter, "TestUser");
            
            % 获取会话节点
            testCase.sessionNode = testCase.protocolNode.getChild(testCase.testSessionName);
            testCase.verifyNotEmpty(testCase.sessionNode);
            
            % 验证会话修改
            testCase.verifyEqual(testCase.sessionNode.sessionInfo.desc, "修改后的会话描述");
            testCase.verifyEqual(testCase.sessionNode.sessionInfo.type, "MEG_RECORDING");
            
            fprintf("✓ 修改后保存并重新加载测试通过\n");
        end
        
        function testDataLoadingAndUnloading(testCase)
            % 测试数据加载和卸载功能
            
            fprintf("\n=== 测试数据加载和卸载功能 ===\n");
            
            % 创建基本项目结构
            testCase.createBasicProjectStructure();
            
            % 1. 初始状态验证
            fprintf("1. 初始状态验证\n");
            testCase.verifyFalse(testCase.protocolNode.isLoaded, "协议初始状态应为未加载");
            testCase.verifyFalse(testCase.protocolNode.chanlocsNode.isLoaded, "通道节点初始状态应为未加载");
            testCase.verifyFalse(testCase.sessionNode.isLoaded, "会话初始状态应为未加载");
            
            % 2. 加载协议数据
            fprintf("2. 加载协议数据\n");
            testCase.protocolNode.load();
            
            testCase.verifyTrue(testCase.protocolNode.isLoaded, "协议应已加载");
            testCase.verifyTrue(testCase.protocolNode.chanlocsNode.isLoaded, "通道节点应已加载");
            testCase.verifyTrue(testCase.protocolNode.cortexNode.isLoaded, "皮层节点应已加载");
            testCase.verifyTrue(testCase.protocolNode.leadfieldNode.isLoaded, "导联场节点应已加载");
            
            % 验证数据已加载
            testCase.verifyNotEmpty(testCase.protocolNode.chanlocsNode.data, "通道数据应已加载");
            testCase.verifyNotEmpty(testCase.protocolNode.cortexNode.data, "皮层数据应已加载");
            testCase.verifyNotEmpty(testCase.protocolNode.leadfieldNode.data, "导联场数据应已加载");
            
            % 3. 加载会话数据
            fprintf("3. 加载会话数据\n");
            testCase.sessionNode.load();
            testCase.verifyTrue(testCase.sessionNode.isLoaded, "会话应已加载");
            
            for i = 1:testCase.sessionNode.childCount
                dataNode = testCase.sessionNode.children(i);
                testCase.verifyTrue(dataNode.isLoaded, sprintf("数据节点 %d 应已加载", i));
                testCase.verifyNotEmpty(dataNode.data, sprintf("数据节点 %d 数据应已加载", i));
            end
            
            % 4. 卸载数据
            fprintf("4. 卸载数据\n");
            testCase.protocolNode.unload();
            testCase.verifyFalse(testCase.protocolNode.isLoaded, "协议应已卸载");
            testCase.verifyFalse(testCase.protocolNode.chanlocsNode.isLoaded, "通道节点应已卸载");
            testCase.verifyFalse(testCase.protocolNode.cortexNode.isLoaded, "皮层节点应已卸载");
            testCase.verifyFalse(testCase.protocolNode.leadfieldNode.isLoaded, "导联场节点应已卸载");
            
            testCase.sessionNode.unload();
            testCase.verifyFalse(testCase.sessionNode.isLoaded, "会话应已卸载");
            
            for i = 1:testCase.sessionNode.childCount
                dataNode = testCase.sessionNode.children(i);
                testCase.verifyFalse(dataNode.isLoaded, sprintf("数据节点 %d 应已卸载", i));
            end
            
            % 5. 重新加载验证
            fprintf("5. 重新加载验证\n");
            testCase.protocolNode.load();
            testCase.sessionNode.load();
            
            testCase.verifyTrue(testCase.protocolNode.isLoaded, "协议应能重新加载");
            testCase.verifyTrue(testCase.sessionNode.isLoaded, "会话应能重新加载");
            
            fprintf("✓ 数据加载和卸载功能测试通过\n");
        end
        
        function testTreeStructureOperations(testCase)
            % 测试树形结构操作
            
            fprintf("\n=== 测试树形结构操作 ===\n");
            
            % 创建基本项目结构
            testCase.createBasicProjectStructure();
            
            % 1. 测试树形结构遍历
            fprintf("1. 测试树形结构遍历\n");
            
            % 打印树形结构
            fprintf("项目树形结构:\n");
            testCase.projectNode.printTree();
            
            % 2. 测试查找子节点
            fprintf("2. 测试查找子节点\n");
            
            % 查找协议节点
            foundProtocol = testCase.projectNode.getChild(testCase.testProtocolName);
            testCase.verifyNotEmpty(foundProtocol);
            testCase.verifyEqual(foundProtocol, testCase.protocolNode);
            
            % 查找会话节点
            foundSession = testCase.protocolNode.getChild(testCase.testSessionName);
            testCase.verifyNotEmpty(foundSession);
            testCase.verifyEqual(foundSession, testCase.sessionNode);
            
            % 3. 测试递归查找
            fprintf("3. 测试递归查找\n");
            
            % 为会话添加第二个数据节点
            testDataPath = fullfile(testCase.testRootPath, "TestData");
            extraDataPath = fullfile(testDataPath, "test_cortex.mat");
            
            extraDataInfo = DataInfo.fromData(extraDataPath, "Name", "extra_data");
            extraDataNode = DataNode.fromInfo(extraDataInfo);
            extraDataNode.path = testCase.sessionNode.path;
            testCase.sessionNode.addChild(extraDataNode);
            
            % 从项目根节点递归查找数据节点
            allDataNodes = testCase.projectNode.findChildren(@(node) isa(node, "DataNode"));
            testCase.verifyEqual(length(allDataNodes), testCase.sessionNode.childCount);
            
            % 4. 测试节点属性
            fprintf("4. 测试节点属性\n");
            
            % 测试根节点属性
            testCase.verifyTrue(testCase.projectNode.isRoot);
            testCase.verifyFalse(testCase.projectNode.isLeaf);
            
            % 测试叶节点属性（数据节点应该是叶节点）
            if testCase.sessionNode.childCount > 0
                dataNode = testCase.sessionNode.children(1);
                testCase.verifyTrue(dataNode.isLeaf);
                testCase.verifyFalse(dataNode.isRoot);
            end
            
            % 5. 测试节点选择
            fprintf("5. 测试节点选择\n");
            
            testCase.protocolNode.select();
            testCase.verifyTrue(testCase.protocolNode.isSelected);
            
            testCase.protocolNode.deselect();
            testCase.verifyFalse(testCase.protocolNode.isSelected);
            
            fprintf("✓ 树形结构操作测试通过\n");
        end
        
        function testFileSystemOperations(testCase)
            % 测试文件系统操作
            
            fprintf("\n=== 测试文件系统操作 ===\n");
            
            % 创建基本项目结构
            testCase.createBasicProjectStructure();
            
            % 1. 验证目录结构
            fprintf("1. 验证目录结构\n");
            testCase.verifyProjectFileStructure();
            
            % 2. 测试从现有数据文件打开
            fprintf("2. 测试从现有数据文件打开\n");
            
            testDataPath = fullfile(testCase.testRootPath, "TestData");
            
            % 测试从现有通道文件打开
            chanlocsDataPath = fullfile(testDataPath, "test_chanlocss.mat");
            chanlocsNode = testCase.protocolNode.openChanlocsFromData(chanlocsDataPath);
            testCase.verifyNotEmpty(chanlocsNode);
            testCase.verifyInstanceOf(chanlocsNode, "ChanlocsNode");
            testCase.verifyEqual(chanlocsNode.name, "test_chanlocss");
            
            % 测试从现有皮层文件打开
            cortexDataPath = fullfile(testDataPath, "test_cortex.mat");
            cortexNode = testCase.protocolNode.openCortexFromData(cortexDataPath);
            testCase.verifyNotEmpty(cortexNode);
            testCase.verifyInstanceOf(cortexNode, "CortexNode");
            testCase.verifyEqual(cortexNode.name, "test_cortex");
            
            % 测试从现有导联场文件打开
            leadfieldDataPath = fullfile(testDataPath, "test_leadfield.mat");
            leadfieldNode = testCase.protocolNode.openLeadfieldFromData(leadfieldDataPath);
            testCase.verifyNotEmpty(leadfieldNode);
            testCase.verifyInstanceOf(leadfieldNode, "LeadfieldNode");
            testCase.verifyEqual(leadfieldNode.name, "test_leadfield");
            
            % 测试从现有数据文件打开
            dataPath = fullfile(testDataPath, "original_data.mat");
            dataNode = testCase.sessionNode.openDataFromData(dataPath);
            testCase.verifyNotEmpty(dataNode);
            testCase.verifyInstanceOf(dataNode, "DataNode");
            testCase.verifyEqual(dataNode.name, "original_data");
            
            % 4. 测试打开现有项目
            fprintf("4. 测试打开现有项目\n");
            
            % 从项目文件打开
            projectFilePath = fullfile(testCase.projectNode.path, testCase.testProjectName + ".mat");
            loadedProject = ProjectNode.openExisting(projectFilePath);
            testCase.verifyNotEmpty(loadedProject);
            testCase.verifyEqual(loadedProject.name, testCase.testProjectName);
            
            % 从项目目录打开
            loadedProject2 = ProjectNode.openExisting(testCase.projectNode.path);
            testCase.verifyNotEmpty(loadedProject2);
            testCase.verifyEqual(loadedProject2.name, testCase.testProjectName);
            
            % 清理
            if isvalid(loadedProject), delete(loadedProject); end
            if isvalid(loadedProject2), delete(loadedProject2); end
            
            fprintf("✓ 文件系统操作测试通过\n");
        end
        
%         function testErrorHandling(testCase)
%             % 测试错误处理
%             
%             fprintf("\n=== 测试错误处理 ===\n");
%             
%             % 1. 测试无效路径
%             fprintf("1. 测试无效路径\n");
%             
%             % 测试打开不存在的项目
%             invalidPath = fullfile(testCase.testRootPath, "NonExistentProject");
%             testCase.verifyError(@() ProjectNode.openExisting(invalidPath), ...
%                 "SEAL:ProjectNode:LoadFailed");
%             
%             % 测试打开不存在的协议
%             testCase.verifyError(@() ProtocolNode.openExisting(invalidPath), ...
%                 "SEAL:ProtocolNode:OpenFailed");
%             
%             % 测试打开不存在的会话
%             testCase.verifyError(@() SessionNode.openExisting(invalidPath), ...
%                 "SEAL:SessionNode:OpenFailed");
%             
%             % 2. 测试无效文件格式
%             fprintf("2. 测试无效文件格式\n");
%             
%             % 创建一个无效的MAT文件
%             invalidFile = fullfile(testCase.testRootPath, "invalid.mat");
%             invalidData = "not a valid structure";
%             save(invalidFile, "invalidData");
%             
%             % 测试从无效文件创建DataInfo
%             testCase.verifyError(@() DataInfo.fromData(invalidFile), ...
%                 "MATLAB:InputParser:ArgumentFailedValidation");
%             
%             % 3. 测试无效节点操作
%             fprintf("3. 测试无效节点操作\n");
%             
%             % 创建基本项目结构
%             testCase.createBasicProjectStructure();
%             
%             % 测试添加错误类型的子节点到协议
%             wrongNode = ChanlocsNode(); % ChanlocsNode不能直接作为ProtocolNode的子节点
%             testCase.verifyError(@() testCase.protocolNode.addChild(wrongNode), ...
%                 "SEAL:ProtocolNode:InvalidChildType");
%             
%             % 测试添加错误类型的子节点到项目
%             testCase.verifyError(@() testCase.projectNode.addChild(wrongNode), ...
%                 "SEAL:ProjectNode:InvalidChildType");
%             
%             % 4. 测试保存到无效路径
%             fprintf("4. 测试保存到无效路径\n");
%             
%             % 测试保存到空路径
%             originalPath = testCase.projectNode.path;
%             testCase.projectNode.path = "";
%             testCase.verifyError(@() testCase.projectNode.save(), ...
%                 "SEAL:ProjectInfo:InvalidPath");
%             
%             % 恢复路径
%             testCase.projectNode.path = originalPath;
%             
%             % 5. 测试加载不存在的数据文件
%             fprintf("5. 测试加载不存在的数据文件\n");
%             
%             % 创建一个无效的DataInfo
%             invalidDataInfo = DataInfo("InvalidData", "non_existent_file.mat", "", "", struct());
%             invalidDataNode = DataNode.fromInfo(invalidDataInfo);
%             
%             % 加载应该失败，但不应崩溃
%             try
%                 invalidDataNode.load();
%                 % 如果到达这里，意味着没有抛出错误
%                 % 根据实际情况，这里可能应该验证某些条件
%             catch ME
%                 fprintf("    预期错误: %s\n", ME.message);
%             end
%             
%             fprintf("✓ 错误处理测试通过\n");
%         end
    end
    
    methods (Access = private)
        function createBasicProjectStructure(testCase)
            % 创建基本的项目结构（项目->协议->会话）
            
            % 创建项目
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testRootPath, ...
                "Description", "基本测试项目");
            
            % 创建协议
            testCase.protocolNode = testCase.projectNode.createNewProtocol(...
                testCase.testProtocolName, ...
                "Type", "EEG_STUDY", ...
                "Description", "基本测试协议");
            
            % 为协议添加单例节点
            testDataPath = fullfile(testCase.testRootPath, "TestData");
            
            chanlocsDataPath = fullfile(testDataPath, "test_chanlocss.mat");
            testCase.protocolNode.openChanlocsFromData(chanlocsDataPath);
            
            cortexDataPath = fullfile(testDataPath, "test_cortex.mat");
            testCase.protocolNode.openCortexFromData(cortexDataPath);
            
            leadfieldDataPath = fullfile(testDataPath, "test_leadfield.mat");
            testCase.protocolNode.openLeadfieldFromData(leadfieldDataPath);
            
            % 创建会话
            testCase.sessionNode = SessionNode.createNew(...
                testCase.testSessionName, ...
                fullfile(testCase.protocolNode.path, "Sessions"), ...
                "Type", "EEG_RECORDING", ...
                "Description", "基本测试会话");

            testCase.protocolNode.addChild(testCase.sessionNode);
            
            % 为会话添加数据节点
            for i = 1:length(testCase.testDataNames)
                dataName = testCase.testDataNames(i);
                dataPath = fullfile(testDataPath, dataName + ".mat");
                
                dataInfo = DataInfo.fromData(dataPath, "Name", dataName);
                dataNode = DataNode.fromInfo(dataInfo);
                dataNode.path = testCase.sessionNode.path;
                testCase.sessionNode.addChild(dataNode);
            end
            
            % 保存项目
            testCase.projectNode.save();
        end
        
        function verifyProjectFileStructure(testCase)
            % 验证项目文件结构
            
            fprintf("验证项目文件结构:\n");
            
            % 项目目录结构
            projectDir = testCase.projectNode.path;
            fprintf("  项目目录: %s\n", projectDir);
            testCase.verifyTrue(isfolder(projectDir), "项目目录应该存在");
            
            % 项目信息文件
            projectInfoFile = fullfile(projectDir, testCase.testProjectName + ".mat");
            fprintf("    项目信息文件: %s\n", projectInfoFile);
            testCase.verifyTrue(isfile(projectInfoFile), "项目信息文件应该存在");
            
            % Protocols目录
            protocolsDir = fullfile(projectDir, "Protocols");
            fprintf("    Protocols目录: %s\n", protocolsDir);
            testCase.verifyTrue(isfolder(protocolsDir), "Protocols目录应该存在");
            
            % 协议目录
            protocolDir = fullfile(protocolsDir, testCase.testProtocolName);
            fprintf("      协议目录: %s\n", protocolDir);
            testCase.verifyTrue(isfolder(protocolDir), "协议目录应该存在");
            
            % 协议信息文件
            protocolInfoFile = fullfile(protocolDir, testCase.testProtocolName + ".mat");
            fprintf("        协议信息文件: %s\n", protocolInfoFile);
            testCase.verifyTrue(isfile(protocolInfoFile), "协议信息文件应该存在");
            
            % 单例节点信息文件
            chanlocsInfoFile = fullfile(protocolDir, "chanlocs_info.mat");
            fprintf("        通道信息文件: %s\n", chanlocsInfoFile);
            testCase.verifyTrue(isfile(chanlocsInfoFile), "通道信息文件应该存在");
            
            cortexInfoFile = fullfile(protocolDir, "cortex_info.mat");
            fprintf("        皮层信息文件: %s\n", cortexInfoFile);
            testCase.verifyTrue(isfile(cortexInfoFile), "皮层信息文件应该存在");
            
            leadfieldInfoFile = fullfile(protocolDir, "leadfield_info.mat");
            fprintf("        导联场信息文件: %s\n", leadfieldInfoFile);
            testCase.verifyTrue(isfile(leadfieldInfoFile), "导联场信息文件应该存在");
            
            % Sessions目录
            sessionsDir = fullfile(protocolDir, "Sessions");
            fprintf("        Sessions目录: %s\n", sessionsDir);
            testCase.verifyTrue(isfolder(sessionsDir), "Sessions目录应该存在");
            
            % 会话目录
            sessionDir = fullfile(sessionsDir, testCase.testSessionName);
            fprintf("          会话目录: %s\n", sessionDir);
            testCase.verifyTrue(isfolder(sessionDir), "会话目录应该存在");
            
            % 会话信息文件
            sessionInfoFile = fullfile(sessionDir, "session_info.mat");
            fprintf("            会话信息文件: %s\n", sessionInfoFile);
            testCase.verifyTrue(isfile(sessionInfoFile), "会话信息文件应该存在");
            
            % 数据文件
            fprintf("            数据文件:\n");
            for i = 1:length(testCase.testDataNames)
                dataFile = fullfile(sessionDir, testCase.testDataNames(i) + ".mat");
                fprintf("              %s\n", dataFile);
                testCase.verifyTrue(isfile(dataFile), sprintf("数据文件 %s 应该存在", testCase.testDataNames(i)));
            end
            
            fprintf("项目文件结构验证通过\n");
        end
    end
    
    methods (Static)
        function runAllTests()
            % 运行所有集成测试
            
            fprintf("开始运行SEAL框架集成测试...\n");
            fprintf("==================================\n");
            
            % 创建测试套件
            suite = matlab.unittest.TestSuite.fromClass(?TestSEALIntegration);
            
            % 运行测试
            results = run(suite);
            
            % 显示结果
            fprintf("\n=== 集成测试结果汇总 ===\n");
            fprintf("运行测试总数: %d\n", numel(results));
            fprintf("通过: %d\n", nnz([results.Passed]));
            fprintf("失败: %d\n", nnz([results.Failed]));
            fprintf("跳过: %d\n", nnz([results.Incomplete]));
            
            if all([results.Passed])
                fprintf("\n✅ 所有集成测试通过！\n");
            else
                fprintf("\n❌ 存在测试失败\n");
            end
            
            fprintf("\n集成测试完成！\n");
        end
        
        function runSelectedTests(testNames)
            % 运行选定的测试
            
            fprintf("开始运行选定的SEAL框架集成测试...\n");
            fprintf("==================================\n");
            
            % 创建测试套件
            suite = matlab.unittest.TestSuite.fromClass(?TestSEALIntegration);
            
            % 筛选测试
            if nargin > 0
                selectedSuite = [];
                for i = 1:length(testNames)
                    testName = testNames{i};
                    if startsWith(testName, "test")
                        testSuite = suite.selectIf("Name", testName);
                        selectedSuite = [selectedSuite, testSuite]; %#ok<AGROW>
                    else
                        % 如果不是以test开头，添加test前缀
                        testSuite = suite.selectIf("Name", "test" + testName);
                        selectedSuite = [selectedSuite, testSuite]; %#ok<AGROW>
                    end
                end
                suite = selectedSuite;
            end
            
            % 运行测试
            results = run(suite);
            
            % 显示结果
            fprintf("\n=== 选定测试结果汇总 ===\n");
            fprintf("运行测试总数: %d\n", numel(results));
            fprintf("通过: %d\n", nnz([results.Passed]));
            fprintf("失败: %d\n", nnz([results.Failed]));
            
            if all([results.Passed])
                fprintf("\n✅ 所有选定测试通过！\n");
            else
                fprintf("\n❌ 存在测试失败\n");
            end
            
            fprintf("\n选定测试完成！\n");
        end
        
        function testPerformance()
            % 性能测试
            
            fprintf("开始SEAL框架性能测试...\n");
            fprintf("=========================\n");
            
            % 创建测试实例
            testCase = TestSEALIntegration();
            testCase.setupClass();
            
            % 测试创建大型项目
            fprintf("\n1. 测试创建大型项目\n");
            
            largeProjectName = "LargePerformanceTest";
            largeProjectPath = fullfile(testCase.testRootPath, "PerformanceTests");
            
            if isfolder(largeProjectPath)
                rmdir(largeProjectPath, "s");
            end
            
            tic;
            largeProject = ProjectNode.createNew(largeProjectName, largeProjectPath);
            
            % 创建多个协议
            nProtocols = 5;
            for i = 1:nProtocols
                protocolName = sprintf("Protocol_%d", i);
                protocol = largeProject.createNewProtocol(protocolName, ...
                    "Type", "EEG_STUDY", "Description", sprintf("性能测试协议 %d", i));
                
                % 添加单例节点
                testDataPath = fullfile(testCase.testRootPath, "TestData");
                protocol.openChanlocsFromData(fullfile(testDataPath, "test_chanlocss.mat"));
                protocol.openCortexFromData(fullfile(testDataPath, "test_cortex.mat"));
                
                % 创建多个会话
                nSessions = 3;
                for j = 1:nSessions
                    sessionName = sprintf("Session_%d_%d", i, j);
                    session = SessionNode.createNew(sessionName, ...
                        fullfile(protocol.path, "Sessions"), ...
                        "Type", "EEG_RECORDING", "Description", sprintf("性能测试会话 %d-%d", i, j));
                    
                    % 添加多个数据节点
                    nDataNodes = 4;
                    for k = 1:nDataNodes
                        dataName = sprintf("Data_%d_%d_%d", i, j, k);
                        dataPath = fullfile(testDataPath, "original_data.mat");
                        dataInfo = DataInfo.fromData(dataPath, "Name", dataName);
                        dataNode = DataNode.fromInfo(dataInfo);
                        dataNode.path = session.path;
                        session.addChild(dataNode);
                    end
                    
                    protocol.addChild(session);
                end
            end
            
            % 保存大型项目
            largeProject.save();
            createTime = toc;
            
            fprintf("  创建大型项目耗时: %.2f 秒\n", createTime);
            fprintf("  创建了 %d 个协议\n", nProtocols);
            fprintf("  每个协议包含约 %d 个会话\n", nSessions);
            fprintf("  每个会话包含约 %d 个数据节点\n", nDataNodes);
            
            % 测试加载性能
            fprintf("\n2. 测试加载性能\n");
            
            tic;
            loadedProject = ProjectNode.openExisting(largeProject.path);
            loadTime = toc;
            
            fprintf("  加载大型项目耗时: %.2f 秒\n", loadTime);
            
            % 测试数据加载性能
            fprintf("\n3. 测试数据加载性能\n");
            
            tic;
            loadedProject.load();
            dataLoadTime = toc;
            
            fprintf("  加载项目所有数据耗时: %.2f 秒\n", dataLoadTime);
            
            % 测试数据卸载性能
            fprintf("\n4. 测试数据卸载性能\n");
            
            tic;
            loadedProject.unload();
            dataUnloadTime = toc;
            
            fprintf("  卸载项目所有数据耗时: %.2f 秒\n", dataUnloadTime);
            
            % 清理
            delete(largeProject);
            delete(loadedProject);
            testCase.teardownClass();
            
            fprintf("\n性能测试完成！\n");
        end
    end
end
classdef TestProtocolNode < matlab.unittest.TestCase
    %TESTPROTOCOLNODE 测试ProtocolNode类
    
    properties
        testProjectName = "TestProject"
        testProtocolName = "TestProtocol"
        testProjectPath
        projectNode
        protocolNode
    end
    
    methods (TestMethodSetup)
        function setup(testCase)
            % 设置测试环境
            testCase.testProjectPath = fullfile(tempdir, "SEAL", 'SEAL_Test_Projects');

            % 确保测试目录存在
            if ~isfolder(testCase.testProjectPath)
                mkdir(testCase.testProjectPath);
            end
            
            % 创建测试项目
            testCase.projectNode = ProjectNode.createNewProject(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
        end
    end
    
    methods (TestMethodTeardown)
        function teardown(testCase)
            % 清理测试环境
            
            % 删除项目节点
            if ~isempty(testCase.projectNode) && isvalid(testCase.projectNode)
                delete(testCase.projectNode);
            end
            
            % 使用更彻底的清理方法
            testCase.cleanupDirectory(testCase.testProjectPath);
        end
    end

    methods (Access = private)
        function cleanupDirectory(testCase, dirPath)
            %CLEANUPDIRECTORY 清理指定目录及其所有内容
            if ~isfolder(dirPath)
                return;
            end
            
            try
                % 删除目录及其所有内容
                rmdir(dirPath, 's');
            catch ME
                % 如果目录删除失败，尝试逐个删除子项
                testCase.cleanupDirectoryContents(dirPath);
            end
        end
    end
    
    methods (Test)
        function testCreateNewProtocol(testCase)
            % 测试创建新协议
            
            fprintf('\n=== 测试创建新协议 ===\n');
            
            % 创建新协议
            testCase.protocolNode = ProtocolNode.createNewProtocol(...
                testCase.testProtocolName, ...
                testCase.projectNode, ...
                'Type', 'EEG', ...
                'Description', '这是一个测试协议');
            
            % 验证协议属性
            testCase.verifyEqual(testCase.protocolNode.name, string(testCase.testProtocolName));
            testCase.verifyEqual(testCase.protocolNode.protocolType, "EEG");
            testCase.verifyEqual(testCase.protocolNode.desc, "这是一个测试协议");
            testCase.verifyEqual(testCase.protocolNode.sessionCount, 0);
            
            % 验证协议文件是否存在
            protocolFile = fullfile(testCase.projectNode.path, 'Protocols', ...
                testCase.testProtocolName, strcat(testCase.testProtocolName, '.mat'));
            testCase.verifyTrue(isfile(protocolFile), '协议文件应该存在');
            
            % 验证项目中的协议数量
            testCase.verifyEqual(testCase.projectNode.protocolCount, 1);
            
            fprintf('✓ 创建新协议测试通过\n');
        end
        
        function testLoadProtocol(testCase)
            % 测试加载协议
            
            fprintf('\n=== 测试加载协议 ===\n');
            
            % 先创建协议
            testCase.protocolNode = ProtocolNode.createNewProtocol(...
                testCase.testProtocolName, ...
                testCase.projectNode);
            
            % 获取协议路径
            protocolPath = fullfile(testCase.projectNode.path, 'Protocols', testCase.testProtocolName);
            
            % 然后加载协议
            loadedProtocol = ProtocolNode.openProtocol(protocolPath, testCase.projectNode);
            
            % 验证加载的协议属性
            testCase.verifyEqual(loadedProtocol.name, string(testCase.testProtocolName));
            testCase.verifyTrue(loadedProtocol.isLoaded);
            
            % 验证父子关系
            testCase.verifyEqual(loadedProtocol.parent, testCase.projectNode);
            
            % 清理
            delete(loadedProtocol);
            
            fprintf('✓ 加载协议测试通过\n');
        end
        
        function testSaveProtocol(testCase)
            % 测试保存协议
            
            fprintf('\n=== 测试保存协议 ===\n');
            
            % 创建协议
            testCase.protocolNode = ProtocolNode.createNewProtocol(...
                testCase.testProtocolName, ...
                testCase.projectNode);
            
            % 修改协议信息
            testCase.protocolNode.protocolInfo.desc = "修改后的协议描述";
            testCase.protocolNode.protocolInfo.type = "MEG";
            
            % 保存协议
            testCase.protocolNode.save();
            
            % 重新加载验证修改是否保存
            protocolPath = fullfile(testCase.projectNode.path, 'Protocols', testCase.testProtocolName);
            reloadedProtocol = ProtocolNode.openProtocol(protocolPath, testCase.projectNode);
            
            testCase.verifyEqual(reloadedProtocol.protocolInfo.desc, "修改后的协议描述");
            testCase.verifyEqual(reloadedProtocol.protocolInfo.type, "MEG");
            
            % 清理
            delete(reloadedProtocol);
            
            fprintf('✓ 保存协议测试通过\n');
        end
        
        function testUnloadProtocol(testCase)
            % 测试卸载协议
            
            fprintf('\n=== 测试卸载协议 ===\n');
            
            % 创建并加载协议
            testCase.protocolNode = ProtocolNode.createNewProtocol(...
                testCase.testProtocolName, ...
                testCase.projectNode);

            % 验证协议已加载
            testCase.verifyTrue(testCase.protocolNode.isLoaded);
            
            % 卸载协议
            testCase.protocolNode.unload();
            
            % 验证协议已卸载
            testCase.verifyFalse(testCase.protocolNode.isLoaded);
            
            fprintf('✓ 卸载协议测试通过\n');
        end
        
        function testProtocolProperties(testCase)
            % 测试协议属性
            
            fprintf('\n=== 测试协议属性 ===\n');
            
            % 创建协议
            testCase.protocolNode = ProtocolNode.createNewProtocol(...
                testCase.testProtocolName, ...
                testCase.projectNode, ...
                'Type', 'EEG', ...
                'Description', '属性测试协议');
            
            % 测试依赖属性
            testCase.verifyEqual(testCase.protocolNode.name, testCase.testProtocolName);
            testCase.verifyEqual(testCase.protocolNode.protocolType, "EEG");
            testCase.verifyEqual(testCase.protocolNode.desc, "属性测试协议");
            testCase.verifyEqual(testCase.protocolNode.sessionCount, 0);
            testCase.verifyEqual(testCase.protocolNode.type, "ProtocolNode");
            
            % 测试协议信息
            testCase.verifyEqual(testCase.protocolNode.protocolInfo.name, testCase.testProtocolName);
            testCase.verifyEqual(testCase.protocolNode.protocolInfo.type, "EEG");
            testCase.verifyEqual(testCase.protocolNode.protocolInfo.desc, "属性测试协议");
            
            fprintf('✓ 协议属性测试通过\n');
        end
        
        function testProjectProtocolRelationship(testCase)
            % 测试项目与协议的关系
            
            fprintf('\n=== 测试项目与协议的关系 ===\n');
            
            % 创建协议
            testCase.protocolNode = ProtocolNode.createNewProtocol(...
                testCase.testProtocolName, ...
                testCase.projectNode);
            
            % 验证父子关系
            testCase.verifyEqual(testCase.protocolNode.parent, testCase.projectNode);
            testCase.verifyEqual(testCase.projectNode.protocolCount, 1);
            
            % 验证项目中的协议列表
            protocols = testCase.projectNode.children;
            testCase.verifyEqual(length(protocols), 1);
            testCase.verifyEqual(protocols(1).name, testCase.testProtocolName);
            
            % 验证协议路径
            expectedPath = fullfile(testCase.projectNode.path, 'Protocols', testCase.testProtocolName);
            testCase.verifyEqual(testCase.protocolNode.path, expectedPath);
            
            fprintf('✓ 项目与协议关系测试通过\n');
        end
        
        function testMultipleProtocols(testCase)
            % 测试多个协议
            
            fprintf('\n=== 测试多个协议 ===\n');
            
            % 创建多个协议
            protocol1 = ProtocolNode.createNewProtocol("Protocol1", testCase.projectNode);
            protocol2 = ProtocolNode.createNewProtocol("Protocol2", testCase.projectNode);
            protocol3 = ProtocolNode.createNewProtocol("Protocol3", testCase.projectNode);
            
            % 验证项目中的协议数量
            testCase.verifyEqual(testCase.projectNode.protocolCount, 3);
            testCase.verifyEqual(testCase.projectNode.childCount, 3);
            
            % 验证每个协议都能正确访问
            protocols = testCase.projectNode.children;
            protocolNames = [protocols.name];
            testCase.verifyTrue(ismember("Protocol1", protocolNames));
            testCase.verifyTrue(ismember("Protocol2", protocolNames));
            testCase.verifyTrue(ismember("Protocol3", protocolNames));
            
            % 清理
            delete(protocol1);
            delete(protocol2);
            delete(protocol3);
            
            fprintf('✓ 多个协议测试通过\n');
        end
        
        function testProtocolDirectoryStructure(testCase)
            % 测试协议目录结构
            
            fprintf('\n=== 测试协议目录结构 ===\n');
            
            % 创建协议
            testCase.protocolNode = ProtocolNode.createNewProtocol(...
                testCase.testProtocolName, ...
                testCase.projectNode);
            
            % 验证目录结构
            protocolDir = fullfile(testCase.projectNode.path, 'Protocols', testCase.testProtocolName);
            sessionsDir = fullfile(protocolDir, 'Sessions');
            
            testCase.verifyTrue(isfolder(protocolDir), '协议目录应该存在');
            testCase.verifyTrue(isfolder(sessionsDir), '会话目录应该存在');
            
            % 验证协议文件
            protocolFile = fullfile(protocolDir, strcat(testCase.testProtocolName, '.mat'));
            testCase.verifyTrue(isfile(protocolFile), '协议文件应该存在');
            
            fprintf('✓ 协议目录结构测试通过\n');
        end
    end
    
    methods (Test, TestTags = {'Integration'})
        function testProjectProtocolIntegration(testCase)
            % 集成测试：项目与协议的完整工作流程
            
            fprintf('\n=== 集成测试：项目与协议的完整工作流程 ===\n');
            
            % 1. 创建项目
            project = ProjectNode.createNewProject(...
                "IntegrationTestProject", ...
                testCase.testProjectPath, ...
                'desc', '集成测试项目');
            
            % 2. 创建多个协议
            protocolEEG = ProtocolNode.createNewProtocol(...
                "EEG_Protocol", project, ...
                'Type', 'EEG', ...
                'Description', 'EEG数据协议');
            
            protocolMEG = ProtocolNode.createNewProtocol(...
                "MEG_Protocol", project, ...
                'Type', 'MEG', ...
                'Description', 'MEG数据协议');
            
            % 3. 验证项目状态
            testCase.verifyEqual(project.protocolCount, 2);
            testCase.verifyEqual(project.name, "IntegrationTestProject");
            
            % 4. 验证协议状态
            testCase.verifyEqual(protocolEEG.protocolType, "EEG");
            testCase.verifyEqual(protocolMEG.protocolType, "MEG");
            testCase.verifyEqual(protocolEEG.parent, project);
            testCase.verifyEqual(protocolMEG.parent, project);
            
            % 5. 保存项目
            project.save();
            
            % 6. 重新加载项目
            projectPath = fullfile(testCase.testProjectPath, "IntegrationTestProject");
            reloadedProject = ProjectNode.openProject(projectPath);
            
            % 7. 验证重新加载的项目
            testCase.verifyEqual(reloadedProject.protocolCount, 2);
            testCase.verifyTrue(reloadedProject.isLoaded);
            
            % 8. 验证重新加载的协议
            protocols = reloadedProject.children;
            testCase.verifyEqual(length(protocols), 2);
            
            eegProtocol = protocols(1);
            megProtocol = protocols(2);
            
            testCase.verifyTrue(ismember("EEG_Protocol", [eegProtocol.name, megProtocol.name]));
            testCase.verifyTrue(ismember("MEG_Protocol", [eegProtocol.name, megProtocol.name]));
            
            % 清理
            delete(protocolEEG);
            delete(protocolMEG);
            delete(project);
            delete(reloadedProject);
            
            fprintf('✓ 项目与协议集成测试通过\n');
        end
    end
    
    methods (Static)
        function runAllTests()
            % 运行所有测试
            
            fprintf('开始运行ProtocolNode测试套件...\n');
            
            % 创建测试套件
            suite = matlab.unittest.TestSuite.fromClass(?TestProtocolNode);
            
            % 运行测试
            results = run(suite);
            
            % 显示结果
            fprintf('\n=== 测试结果 ===\n');
            fprintf('运行测试: %d\n', numel(results));
            fprintf('通过: %d\n', nnz([results.Passed]));
            fprintf('失败: %d\n', nnz([results.Failed]));
            fprintf('跳过: %d\n', nnz([results.Incomplete]));
            
            if all([results.Passed])
                fprintf('\n所有测试通过！\n');
            else
                fprintf('\n有测试失败\n');
            end
        end
        
        function runIntegrationTests()
            % 运行集成测试
            
            fprintf('开始运行ProtocolNode集成测试...\n');
            
            % 创建集成测试套件
            suite = matlab.unittest.TestSuite.fromClass(?TestProtocolNode, 'Tag', 'Integration');
            
            % 运行测试
            results = run(suite);
            
            % 显示结果
            fprintf('\n=== 集成测试结果 ===\n');
            fprintf('运行集成测试: %d\n', numel(results));
            fprintf('通过: %d\n', nnz([results.Passed]));
            fprintf('失败: %d\n', nnz([results.Failed]));
            
            if all([results.Passed])
                fprintf('\n所有集成测试通过！\n');
            else
                fprintf('\n有集成测试失败\n');
            end
        end
    end
end
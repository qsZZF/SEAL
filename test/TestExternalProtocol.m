classdef TestExternalProtocol < matlab.unittest.TestCase
    %TESTEXTERNALPROTOCOL 测试打开外部Protocol并保存的功能
    
    properties
        testProjectName = "TestProject"
        testExternalProtocolName = "ExternalProtocol"
        testProjectPath
        externalProtocolPath
        projectNode
    end
    
    methods (TestMethodSetup)
        function setup(testCase)
            % 设置测试环境
            testCase.testProjectPath = fullfile(tempdir, "SEAL", 'SEAL_Test_Projects');
            testCase.externalProtocolPath = fullfile(tempdir, "SEAL", 'SEAL_External_Protocols');
            
            % 确保测试目录存在
            if ~isfolder(testCase.testProjectPath)
                mkdir(testCase.testProjectPath);
            end
            if ~isfolder(testCase.externalProtocolPath)
                mkdir(testCase.externalProtocolPath);
            end
            
            fprintf('测试目录: %s\n', testCase.testProjectPath);
            fprintf('外部协议目录: %s\n', testCase.externalProtocolPath);
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
            testCase.cleanupDirectory(testCase.externalProtocolPath);
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
        function testOpenExternalProtocolFromFolder(testCase)
            % 测试从文件夹打开外部协议
            
            fprintf('\n=== 测试从文件夹打开外部协议 ===\n');
            
            % 1. 创建项目
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 2. 创建外部协议文件夹结构
            externalProtocolDir = fullfile(testCase.externalProtocolPath, testCase.testExternalProtocolName);
            if ~isfolder(externalProtocolDir)
                mkdir(externalProtocolDir);
            end
            
            % 创建外部协议文件
            externalProtocolInfo = ProtocolInfo.createNew(...
                testCase.testExternalProtocolName, ...
                'EEG', ...
                '这是一个外部协议');
            externalProtocolInfo.save(externalProtocolDir);
            
            % 3. 打开外部协议
            externalProtocol = testCase.projectNode.openProtocol(externalProtocolDir);
            
            % 4. 验证协议属性
            testCase.verifyEqual(externalProtocol.name, testCase.testExternalProtocolName);
            testCase.verifyEqual(externalProtocol.protocolType, "EEG");
            testCase.verifyTrue(externalProtocol.isLoaded);
            
            % 5. 验证协议已添加到项目
            testCase.verifyEqual(testCase.projectNode.protocolCount, 1);
            testCase.verifyEqual(testCase.projectNode.children(1).name, testCase.testExternalProtocolName);
            
            fprintf('✓ 从文件夹打开外部协议测试通过\n');
        end
        
        function testOpenExternalProtocolFromMatFile(testCase)
            % 测试从.mat文件打开外部协议
            
            fprintf('\n=== 测试从.mat文件打开外部协议 ===\n');
            
            % 1. 创建项目
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 2. 创建外部协议文件夹和文件
            externalProtocolDir = fullfile(testCase.externalProtocolPath, testCase.testExternalProtocolName);
            if ~isfolder(externalProtocolDir)
                mkdir(externalProtocolDir);
            end
            
            externalProtocolFile = fullfile(externalProtocolDir, strcat(testCase.testExternalProtocolName, '.mat'));
            externalProtocolInfo = ProtocolInfo.createNew(...
                testCase.testExternalProtocolName, ...
                'MEG', ...
                '从.mat文件打开的外部协议');
            externalProtocolInfo.save(externalProtocolDir);
            
            % 3. 直接从.mat文件打开外部协议
            externalProtocol = testCase.projectNode.openProtocol(externalProtocolFile);
            
            % 4. 验证协议属性
            testCase.verifyEqual(externalProtocol.name, testCase.testExternalProtocolName);
            testCase.verifyEqual(externalProtocol.protocolType, "MEG");
            testCase.verifyTrue(externalProtocol.isLoaded);
            
            % 5. 验证协议已添加到项目
            testCase.verifyEqual(testCase.projectNode.protocolCount, 1);
            
            fprintf('✓ 从.mat文件打开外部协议测试通过\n');
        end
        
        function testSaveExternalProtocolToProject(testCase)
            % 测试将外部协议保存到项目中
            
            fprintf('\n=== 测试将外部协议保存到项目中 ===\n');
            
            % 1. 创建项目
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 2. 创建外部协议
            externalProtocolDir = fullfile(testCase.externalProtocolPath, testCase.testExternalProtocolName);
            if ~isfolder(externalProtocolDir)
                mkdir(externalProtocolDir);
            end
            
            externalProtocolInfo = ProtocolInfo.createNew(...
                testCase.testExternalProtocolName, ...
                'EEG', ...
                '需要保存到项目的外部协议');
            externalProtocolInfo.save(externalProtocolDir);
            
            % 3. 打开外部协议
            externalProtocol = testCase.projectNode.openProtocol(externalProtocolDir);
            
            % 4. 修改协议信息
            externalProtocol.protocolInfo.desc = "修改后的协议描述";
            
            % 5. 保存协议到项目
            externalProtocol.save();
            
            % 6. 验证协议文件已保存到项目目录
            projectProtocolPath = fullfile(testCase.projectNode.path, 'Protocols', testCase.testExternalProtocolName);
            projectProtocolFile = fullfile(projectProtocolPath, strcat(testCase.testExternalProtocolName, '.mat'));
            
            testCase.verifyTrue(isfolder(projectProtocolPath), '协议目录应该在项目内创建');
            testCase.verifyTrue(isfile(projectProtocolFile), '协议文件应该在项目内创建');
            
            % 7. 重新加载验证修改是否保存
            reloadedProtocol = testCase.projectNode.openProtocol(projectProtocolPath);
            testCase.verifyEqual(reloadedProtocol.protocolInfo.desc, "修改后的协议描述");
            
            % 清理
            delete(reloadedProtocol);
            
            fprintf('✓ 外部协议保存到项目测试通过\n');
        end
        
        function testExternalProtocolPathConsistency(testCase)
            % 测试外部协议路径一致性
            
            fprintf('\n=== 测试外部协议路径一致性 ===\n');
            
            % 1. 创建项目
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 2. 创建外部协议
            externalProtocolDir = fullfile(testCase.externalProtocolPath, testCase.testExternalProtocolName);
            if ~isfolder(externalProtocolDir)
                mkdir(externalProtocolDir);
            end
            
            externalProtocolInfo = ProtocolInfo.createNew(...
                testCase.testExternalProtocolName, ...
                'EEG', ...
                '路径一致性测试协议');
            externalProtocolInfo.save(externalProtocolDir);
            
            % 3. 打开外部协议
            externalProtocol = testCase.projectNode.openProtocol(externalProtocolDir);
            
            % 4. 验证路径设置正确
            expectedPath = fullfile(testCase.projectNode.path, 'Protocols', testCase.testExternalProtocolName);
            testCase.verifyEqual(externalProtocol.path, expectedPath);
            
            % 5. 验证保存后路径不变
            externalProtocol.save();
            testCase.verifyEqual(externalProtocol.path, expectedPath);
            
            fprintf('✓ 外部协议路径一致性测试通过\n');
        end
        
        function testMultipleExternalProtocols(testCase)
            % 测试打开多个外部协议
            
            fprintf('\n=== 测试打开多个外部协议 ===\n');
            
            % 1. 创建项目
            testCase.projectNode = ProjectNode.createNew(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 2. 创建多个外部协议
            protocolNames = {"ExternalProtocol1", "ExternalProtocol2", "ExternalProtocol3"};
            protocolTypes = {"EEG", "MEG", "fMRI"};
            
            for i = 1:length(protocolNames)
                protocolDir = fullfile(testCase.externalProtocolPath, protocolNames{i});
                if ~isfolder(protocolDir)
                    mkdir(protocolDir);
                end
                
                protocolInfo = ProtocolInfo.createNew(...
                    protocolNames{i}, ...
                    protocolTypes{i}, ...
                    sprintf('外部协议 %d', i));
                protocolInfo.save(protocolDir);
                
                % 打开外部协议
                testCase.projectNode.openProtocol(protocolDir);
            end
            
            % 3. 验证所有协议都已加载
            testCase.verifyEqual(testCase.projectNode.protocolCount, 3);
            
            % 4. 验证协议属性
            protocols = testCase.projectNode.children;
            for i = 1:length(protocols)
                testCase.verifyEqual(protocols(i).name, protocolNames{i});
                testCase.verifyEqual(protocols(i).protocolType, protocolTypes{i});
            end
            
            fprintf('✓ 多个外部协议测试通过\n');
        end
    end
    
    methods (Test, TestTags = {'Integration'})
        function testExternalProtocolFullWorkflow(testCase)
            % 集成测试：外部协议的完整工作流程
            
            fprintf('\n=== 集成测试：外部协议的完整工作流程 ===\n');
            
            % 1. 创建项目和外部协议
            project = ProjectNode.createNew("WorkflowTestProject", testCase.testProjectPath);
            
            externalProtocolDir = fullfile(testCase.externalProtocolPath, "WorkflowExternalProtocol");
            if ~isfolder(externalProtocolDir)
                mkdir(externalProtocolDir);
            end
            
            externalProtocolInfo = ProtocolInfo.createNew(...
                "WorkflowExternalProtocol", ...
                "EEG", ...
                "完整工作流程测试协议");
            externalProtocolInfo.save(externalProtocolDir);
            
            % 2. 打开外部协议
            externalProtocol = project.openProtocol(externalProtocolDir);
            
            % 3. 修改并保存
            externalProtocol.protocolInfo.desc = "工作流程修改后的描述";
            externalProtocol.save();
            
            % 4. 重新加载项目验证持久化
            projectPath = fullfile(testCase.testProjectPath, "WorkflowTestProject");
            reloadedProject = ProjectNode.openExisting(projectPath);
            
            % 5. 验证外部协议已正确集成到项目中
            testCase.verifyEqual(reloadedProject.protocolCount, 1);
            reloadedProtocol = reloadedProject.children(1);
            testCase.verifyEqual(reloadedProtocol.name, "WorkflowExternalProtocol");
            testCase.verifyEqual(reloadedProtocol.protocolInfo.desc, "工作流程修改后的描述");
            
            % 清理
            delete(project);
            delete(reloadedProject);
            
            fprintf('✓ 外部协议完整工作流程测试通过\n');
        end
    end
    
    methods (Static)
        function runAllTests()
            % 运行所有测试
            
            fprintf('开始运行ExternalProtocol测试套件...\n');
            
            % 创建测试套件
            suite = matlab.unittest.TestSuite.fromClass(?TestExternalProtocol);
            
            % 运行测试
            results = run(suite);
            
            % 显示结果
            fprintf('\n=== 测试结果 ===\n');
            fprintf('运行测试: %d\n', numel(results));
            fprintf('通过: %d\n', nnz([results.Passed]));
            fprintf('失败: %d\n', nnz([results.Failed]));
            fprintf('跳过: %d\n', nnz([results.Incomplete]));
            
            if all([results.Passed])
                fprintf('\n所有外部协议测试通过！\n');
            else
                fprintf('\n测试失败\n');
            end
        end
    end
end
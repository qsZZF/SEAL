classdef TestProjectNode < matlab.unittest.TestCase
    %TESTPROJECTNODE 测试ProjectNode类
    
    properties
        testProjectName = "TestProject"
        testProjectPath
        projectNode
    end
    
    methods (TestMethodSetup)
        function setup(testCase)
            % 设置测试环境
            testCase.testProjectPath = fullfile(tempdir, 'SEAL_Test_Projects');

            % 确保测试目录存在
            if ~isfolder(testCase.testProjectPath)
                mkdir(testCase.testProjectPath);
            end
            
            fprintf('测试目录: %s\n', testCase.testProjectPath);
        end
    end
    
    methods (TestMethodTeardown)
        function teardown(testCase)
            % 清理测试环境
            if ~isempty(testCase.projectNode) && isvalid(testCase.projectNode)
                delete(testCase.projectNode);
            end
            
            % 删除测试项目文件
            if isfolder(fullfile(testCase.testProjectPath, testCase.testProjectName))
                rmdir(fullfile(testCase.testProjectPath, testCase.testProjectName), 's');
            end
        end
    end
    
    methods (Test)
        function testCreateNewProject(testCase)
            % 测试创建新项目
            
            fprintf('\n=== 测试创建新项目 ===\n');
            
            % 创建新项目
            testCase.projectNode = ProjectNode.createNewProject(...
                testCase.testProjectName, ...
                testCase.testProjectPath, ...
                'desc', '这是一个测试项目');
            
            % 验证项目属性
            testCase.verifyEqual(testCase.projectNode.name, string(testCase.testProjectName));
            testCase.verifyEqual(testCase.projectNode.protocolCount, 0);
            
            % 验证项目文件是否存在
            projectFile = fullfile(testCase.testProjectPath, testCase.testProjectName, ...
                strcat(testCase.testProjectName, '.mat'));
            testCase.verifyTrue(isfile(projectFile), '项目文件应该存在');
            
            fprintf('✓ 创建新项目测试通过\n');
        end
        
        function testLoadProject(testCase)
            % 测试加载项目
            
            fprintf('\n=== 测试加载项目 ===\n');
            
            % 先创建项目
            testCase.projectNode = ProjectNode.createNewProject(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 然后加载项目
            projectPath = fullfile(testCase.testProjectPath, testCase.testProjectName);
            loadedProject = ProjectNode.openProject(projectPath);
            
            % 验证加载的项目属性
            testCase.verifyEqual(loadedProject.name, string(testCase.testProjectName));
            testCase.verifyTrue(loadedProject.isLoaded);
            
            % 清理
            delete(loadedProject);
            
            fprintf('✓ 加载项目测试通过\n');
        end
        
%         function testAddProtocol(testCase)
%             % 测试添加协议
%             
%             fprintf('\n=== 测试添加协议 ===\n');
%             
%             % 创建项目
%             testCase.projectNode = TestProjectNode.createNewProject(...
%                 testCase.testProjectName, ...
%                 testCase.testProjectPath);
%             
%             % 创建协议节点（需要ProtocolNode类存在）
%             try
%                 % 这里假设ProtocolNode类存在
%                 protocol = ProtocolNode('TestProtocol');
%                 testCase.projectNode.addChild(protocol);
%                 
%                 % 验证协议数量
%                 testCase.verifyEqual(testCase.projectNode.protocolCount, 1);
%                 testCase.verifyEqual(testCase.projectNode.childCount, 1);
%                 
%                 fprintf('✓ 添加协议测试通过\n');
%                 
%             catch ME
%                 if strcmp(ME.identifier, 'MATLAB:undefinedVarOrClass')
%                     fprintf('⚠  ProtocolNode类未定义，跳过协议添加测试\n');
%                 else
%                     rethrow(ME);
%                 end
%             end
%         end
%         
        function testSaveProject(testCase)
            % 测试保存项目
            
            fprintf('\n=== 测试保存项目 ===\n');
            
            % 创建项目
            testCase.projectNode = ProjectNode.createNewProject(...
                testCase.testProjectName, ...
                testCase.testProjectPath);
            
            % 修改项目信息
            testCase.projectNode.projectInfo.desc = "修改后的描述";
            
            % 保存项目
            testCase.projectNode.save();
            
            % 重新加载验证修改是否保存
            projectPath = fullfile(testCase.testProjectPath, testCase.testProjectName);
            reloadedProject = ProjectNode.openProject(projectPath);
            
            testCase.verifyEqual(reloadedProject.projectInfo.desc, "修改后的描述");
            
            % 清理
            delete(reloadedProject);
            
            fprintf('✓ 保存项目测试通过\n');
        end
        
        function testUnloadProject(testCase)
            % 测试卸载项目
            
            fprintf('\n=== 测试卸载项目 ===\n');
            
            % 创建并加载项目
            testCase.projectNode = ProjectNode.createNewProject(...
                testCase.testProjectName, ...
                testCase.testProjectPath);

            % 验证项目已加载
            testCase.verifyTrue(testCase.projectNode.isLoaded);
            
            % 卸载项目
            testCase.projectNode.unload();
            
            % 验证项目已卸载
            testCase.verifyFalse(testCase.projectNode.isLoaded);
            
            fprintf('✓ 卸载项目测试通过\n');
        end
        
        function testProjectProperties(testCase)
            % 测试项目属性
            
            fprintf('\n=== 测试项目属性 ===\n');
            
            % 创建项目
            testCase.projectNode = ProjectNode.createNewProject(...
                testCase.testProjectName, ...
                testCase.testProjectPath, ...
                'desc', '属性测试项目');
            
            % 测试依赖属性
            testCase.verifyEqual(testCase.projectNode.name, testCase.testProjectName);
            testCase.verifyEqual(testCase.projectNode.protocolCount, 0);
            testCase.verifyEqual(testCase.projectNode.type, "ProjectNode");
            
            % 测试项目信息
            testCase.verifyEqual(testCase.projectNode.projectInfo.name, testCase.testProjectName);
            testCase.verifyEqual(testCase.projectNode.projectInfo.desc, "属性测试项目");
            
            fprintf('✓ 项目属性测试通过\n');
        end
    end
    
    methods (Static)
        function runAllTests()
            % 运行所有测试
            
            fprintf('开始运行ProjectNode测试套件...\n');
            
            % 创建测试套件
            suite = matlab.unittest.TestSuite.fromClass(?TestProjectNode);
            
            % 运行测试
            results = run(suite);
            
            % 显示结果
            fprintf('\n=== 测试结果 ===\n');
            fprintf('运行测试: %d\n', numel(results));
            fprintf('通过: %d\n', nnz([results.Passed]));
            fprintf('失败: %d\n', nnz([results.Failed]));
            fprintf('跳过: %d\n', nnz([results.Incomplete]));
            
            if all([results.Passed])
                fprintf('\n 所有测试通过！\n');
            else
                fprintf('\n❌ 有测试失败，请检查代码\n');
            end
        end
    end
end
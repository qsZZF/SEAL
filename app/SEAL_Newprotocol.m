classdef SEAL_Newprotocol < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        CreateNewProtocolPanel         matlab.ui.container.Panel
        NewProtocolNameEditField       matlab.ui.control.EditField
        NewProtocolNameEditFieldLabel  matlab.ui.control.Label
        CreateButton                   matlab.ui.control.Button
        CancelButton                   matlab.ui.control.Button
    end


    properties (Access = private)
        Callingapp
        ProtocolName = ''
        ProtocolType = 'Real Data Analysis' % Default value
        WasCancelled = true % Assume cancelled until 'Create' is pressed
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.Callingapp=mainapp;
            app.Callingapp.inode = app.Callingapp.inode+1;
            
            protocols = fullfile(app.Callingapp.ProjectPath, 'Protocols');
            defaultProtocolName = findAvailableName(protocols, 'Protocol #');
            if isprop(app, 'NewProtocolNameEditField') && isvalid(app.NewProtocolNameEditField)
                app.NewProtocolNameEditField.Value = defaultProtocolName;
            end
            
            
%             protocolName = sprintf('Protocol_%d', app.Callingapp.inode);
            app.NewProtocolNameEditField.Value=defaultProtocolName;
           if ~isempty(mainapp)
                mainFigPos = mainapp.UIFigure.Position;
                projectPanelLeft = 4; % Corresponds to ProjectPanel position
                dialogWidth = 258;
                dialogHeight = 138;
                newX = mainFigPos(1) + projectPanelLeft;
                newY = mainFigPos(2) + (mainFigPos(4) - dialogHeight) / 2;
                app.UIFigure.Position = [newX, newY, dialogWidth, dialogHeight];
            else
                app.UIFigure.Position = [100 100 258 205];
            end
        end

        % Button pushed function: CreateButton
        function CreateButtonPushed(app, event)
%             protocolName = sprintf('Protocol_%d', app.Callingapp.i);
            
            protocolName = app.NewProtocolNameEditField.Value;
            if isempty(strtrim(protocolName))
                errordlg('Protocol name cannot be empty.', 'Invalid Protocol Name');
                return;
            end
            
            newProtocolPath = fullfile(app.Callingapp.ProjectPath, protocolName);

            if exist(newProtocolPath, 'dir')
                errordlg(['A protocol named "', protocolName, '" already exists in the current base directory.'], 'Error: Protocol Exists');
                return;
            end

            protocolNode = app.Callingapp.projectNode.createNewProtocol(protocolName);

            app.Callingapp.flushFileTree()

            % create directory for new protocol
            delete(app);
            
        end

        % Button pushed function: CancelButton
        function CancelButtonPushed(app, event)
            app.Callingapp.inode = app.Callingapp.inode-1;
            
            app.Callingapp.appendOutput('New protocol creation cancelled.');
            delete(app);
%             uiwait(msgbox('New protocol creation cancelled.'));
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 258 138];
            app.UIFigure.Name = 'MATLAB App';

            % Create CreateNewProtocolPanel
            app.CreateNewProtocolPanel = uipanel(app.UIFigure);
            app.CreateNewProtocolPanel.Title = 'Create New Protocol';
            app.CreateNewProtocolPanel.Position = [1 1 258 138];

            % Create CancelButton
            app.CancelButton = uibutton(app.CreateNewProtocolPanel, 'push');
            app.CancelButton.ButtonPushedFcn = createCallbackFcn(app, @CancelButtonPushed, true);
            app.CancelButton.Position = [133 12 92 22];
            app.CancelButton.Text = 'Cancel';

            % Create CreateButton
            app.CreateButton = uibutton(app.CreateNewProtocolPanel, 'push');
            app.CreateButton.ButtonPushedFcn = createCallbackFcn(app, @CreateButtonPushed, true);
            app.CreateButton.Position = [31 12 92 22];
            app.CreateButton.Text = 'Create';

            % Create NewProtocolNameEditFieldLabel
            app.NewProtocolNameEditFieldLabel = uilabel(app.CreateNewProtocolPanel);
            app.NewProtocolNameEditFieldLabel.HorizontalAlignment = 'right';
            app.NewProtocolNameEditFieldLabel.Position = [22 88 112 22];
            app.NewProtocolNameEditFieldLabel.Text = 'New Protocol Name';

            % Create NewProtocolNameEditField
            app.NewProtocolNameEditField = uieditfield(app.CreateNewProtocolPanel, 'text');
            app.NewProtocolNameEditField.Position = [31 58 194 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SEAL_Newprotocol(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
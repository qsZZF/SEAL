<<<<<<< HEAD
<<<<<<< HEAD
function [seal_struct, metadata] = seal_importer(filePath, dataType, varargin)
%SEAL_IMPORTER Core engine for importing data into SEAL format.
%   This function converts data from various formats into the standardized
%   SEAL Toolbox struct format. It is designed to be called directly from

%   scripts for batch processing.
%
%   Usage:
%   seal_struct = seal_importer_engine(filePath, fileFormat, dataType)
%   seal_struct = seal_importer_engine(filePath, 'generic_mat', dataType, mappings)
%
%   Inputs:
%   - filePath:   Full path to the source data file.
%   - fileFormat: String specifying the format (e.g., 'eeglab_set', 'bst_cortex', 'generic_mat').
%   - dataType:   The target SEAL data type ('Data', 'Cortex', 'LeadField').

%
%   Output:
%   - seal_struct: A standardized struct conforming to the SEAL data template.
metadata = struct('ID', [], 'DataType', [], 'DataSize', [], 'Chanlocs', [],...
                'SamplingRate', [], 'Cortex', [], 'LeadField', [], 'HeadModelType', [], 'Orientation', []);
Data_mapping = struct('ID', 'Data', 'data',[], 'Srate', [], 'Time', [], 'Chanlocs', [], 'Event', [], 'Type', 'EEG');
Cortex_mapping = struct('ID', 'Cortex  ', 'Vertices', [], 'Faces', [], 'VertConn', [], 'Atlas', [], 'Structure', [],'VertNormals',[], 'Type', 'Cortex');
LeadField_mapping = struct('ID', 'LeadField', 'Gain', [], 'Orientation', [],'GridLoc', [], 'GridOrient', [],...
    'HeadModelType', 'Surface', 'Type', 'LeadField');
Chanlocs_mapping = struct('ID', 'Chanlocs', 'chanlocs', [],'Type', 'Chanlocs');
default_mappings = struct('Data', Data_mapping, 'Cortex', Cortex_mapping, 'LeadField', LeadField_mapping,'Chanlocs',Chanlocs_mapping);

    % --- Input Validation ---
    p = inputParser;
    addRequired(p, 'filePath', @ischar);
%     addRequired(p, 'fileFormat', @ischar);
    addRequired(p, 'dataType', @ischar);
    addParameter(p, 'mappings', default_mappings.(dataType), @(x) isstruct(x));
    parse(p, filePath, dataType, varargin{:});
    
    mappings = p.Results.mappings;
    sourceData = load(filePath); 
    if numel(sourceData)==1
        sourceData = struct2cell(sourceData);
        sourceData = sourceData{1};
    end
    if ismatrix(sourceData)&&~isstruct(sourceData)
        sourceData = load(filePath);
    end
    seal_struct = map_from_generic_mat(sourceData, mappings , dataType);


    if ~isfield(seal_struct,'Type')
        seal_struct.Type = default_mappings.(dataType).Type;
    end
    switch dataType
        case 'Data'
            metadata.DataType = seal_struct.Type;
            metadata.DataSize = size(seal_struct.data);
            metadata.SamplingRate = seal_struct.srate;
            if isfield(seal_struct,'chanlocs')
                metadata.Chanlocs = 'Yes';
            end
        case 'Cortex'
            metadata.Cortex = 'Yes';
        case 'LeadField'
            metadata.LeadField = 'Yes';
            metadata.HeadModelType = seal_struct.HeadModelType;
            metadata.Orientation = seal_struct.Orientation;
        case 'Chanlocs'
            metadata.Chanlocs = 'Yes';
    end
        
    % --- Main Switch for different file formats ---
%     switch fileFormat
%         case 'eeglab_set'
%             if ~strcmp(dataType, 'Data')
%                 error('SEAL:ImportError', 'EEGLAB .set format is only applicable for EEG Data.');
%             end
%             seal_struct = import_from_eeglab(filePath);
%             
%         case 'bst_cortex'
%              if ~strcmp(dataType, 'Cortex')
%                 error('SEAL:ImportError', 'Brainstorm surface format is only applicable for Cortex data.');
%             end
%             seal_struct = import_from_bst_surface(filePath);
% 
%         case 'generic_mat'
%             sourceData = load(filePath); sourceData = sourceData{1};
%             seal_struct = map_from_generic_mat(sourceData, mappings);
%             
%         otherwise
%             error('SEAL:ImportError', 'Unsupported file format: %s', fileFormat);
%     end
end

% --- Private Helper Functions ---

function EEG_Data = import_from_eeglab(filePath)
    % Converts an EEGLAB .set file to the SEAL EEG_Data structure
    if ~exist('pop_loadset', 'file')
        error('SEAL:DependencyError', 'EEGLAB function pop_loadset not found. Is EEGLAB in your MATLAB path?');
    end
    EEG = pop_loadset('filename', filePath);
    
    EEG_Data = struct();
    EEG_Data.Comment = EEG.setname;
    EEG_Data.F = EEG.data;
    EEG_Data.srate = EEG.srate;
    EEG_Data.time = EEG.times / 1000; % Convert ms to seconds
    EEG_Data.chanlocs = EEG.chanlocs;
    EEG_Data.events = EEG.event;
    EEG_Data.type = 'EEG';
    EEG_Data.ica = struct();
end

function Cortex = import_from_bst_surface(filePath)
    % Converts a Brainstorm surface .mat file to the SEAL Cortex structure
    s = load(filePath);
    Cortex = struct();
    Cortex.Comment = s.Comment;
    Cortex.Vertices = s.Vertices;
    Cortex.Faces = s.Faces;
    if isfield(s, 'VertConn'), Cortex.VertConn = s.VertConn; end
end

function MappedStruct = map_from_generic_mat(sourceData, mappings, dataType)
    % Applies user-defined mappings to a generic struct
    MappedStruct = struct();
    allMapFields = fieldnames(mappings);

    % --- Process all provided mappings ---
    for i = 1:numel(allMapFields)
        sealField = allMapFields{i};
        mappingValue = mappings.(sealField);

        if ischar(mappingValue) && isfield(sourceData,mappingValue)
            % If the mapping is a string, treat it as a variable name (potentially nested)
            try
                MappedStruct.(sealField) = get_nested_field(sourceData, mappingValue);
            catch ME
                warning('SEAL:MappingWarning', 'Could not resolve source variable ''%s'' for field ''%s''. It will be skipped. Error: %s', mappingValue, sealField, ME.message);
            end
        else
            % If the mapping is not a string, assign the literal value
            MappedStruct.(sealField) = mappingValue;
        end
    end
    
    % --- Validate that all REQUIRED fields now exist in the MappedStruct ---
    requiredFields = {};
    switch dataType
        case 'Data', requiredFields = {'data', 'srate'};
        case 'Cortex', requiredFields = {'Vertices', 'Faces'};
        case 'LeadField', requiredFields = {'Gain', 'Orientation'};
    end
    
    for i = 1:numel(requiredFields)
        field = requiredFields{i};
        if ~isfield(MappedStruct, field)
            error('SEAL:MappingError', 'Required field ''%s'' was not successfully mapped and is missing from the final struct.', field);
        end
    end
    
    % --- Add Default Values ---
%     if ~isfield(MappedStruct, 'ID')
%         MappedStruct.ID = [dataType,'_default'];
%     end
end

function value = get_nested_field(s, fieldvalue)

    
    if isfield(s,fieldvalue)
        value = s.(fieldvalue);
    end
%     parts = strsplit(fieldPath, '.');
%     currentValue = s;
%     for k = 1:numel(parts)
%         if isfield(currentValue, parts{k})
%             currentValue = currentValue.(parts{k});
%         else
%             fprintf('Field path part ''%s'' not found.', parts{k});
%         end
%     end
%     value = currentValue;
=======
function [seal_struct, metadata] = seal_importer(filePath, dataType, varargin)
%SEAL_IMPORTER Core engine for importing data into SEAL format.
%   This function converts data from various formats into the standardized
%   SEAL Toolbox struct format. It is designed to be called directly from

%   scripts for batch processing.
%
%   Usage:
%   seal_struct = seal_importer_engine(filePath, fileFormat, dataType)
%   seal_struct = seal_importer_engine(filePath, 'generic_mat', dataType, mappings)
%
%   Inputs:
%   - filePath:   Full path to the source data file.
%   - fileFormat: String specifying the format (e.g., 'eeglab_set', 'bst_cortex', 'generic_mat').
%   - dataType:   The target SEAL data type ('Data', 'Cortex', 'LeadField').

%
%   Output:
%   - seal_struct: A standardized struct conforming to the SEAL data template.
metadata = struct('ID', [], 'DataType', [], 'DataSize', [], 'Chanlocs', [],...
                'SamplingRate', [], 'Cortex', [], 'LeadField', [], 'HeadModelType', [], 'Orientation', []);
Data_mapping = struct('ID', 'Data', 'data',[], 'Srate', [], 'Time', [], 'Chanlocs', [], 'Event', [], 'Type', 'EEG');
Cortex_mapping = struct('ID', 'Cortex  ', 'Vertices', [], 'Faces', [], 'VertConn', [], 'Atlas', [], 'Structure', [],'VertNormals',[], 'Type', 'Cortex');
LeadField_mapping = struct('ID', 'LeadField', 'Gain', [], 'Orientation', [],'GridLoc', [], 'GridOrient', [],...
    'HeadModelType', 'Surface', 'Type', 'LeadField');
Chanlocs_mapping = struct('ID', 'Chanlocs', 'chanlocs', [],'Type', 'Chanlocs');
default_mappings = struct('Data', Data_mapping, 'Cortex', Cortex_mapping, 'LeadField', LeadField_mapping,'Chanlocs',Chanlocs_mapping);

    % --- Input Validation ---
    p = inputParser;
    addRequired(p, 'filePath', @ischar);
%     addRequired(p, 'fileFormat', @ischar);
    addRequired(p, 'dataType', @ischar);
    addParameter(p, 'mappings', default_mappings.(dataType), @(x) isstruct(x));
    parse(p, filePath, dataType, varargin{:});
    
    mappings = p.Results.mappings;
    sourceData = load(filePath); 
    if numel(sourceData)==1
        sourceData = struct2cell(sourceData);
        sourceData = sourceData{1};
    end
    if ismatrix(sourceData)&&~isstruct(sourceData)
        sourceData = load(filePath);
    end
    seal_struct = map_from_generic_mat(sourceData, mappings , dataType);


    if ~isfield(seal_struct,'Type')
        seal_struct.Type = default_mappings.(dataType).Type;
    end
    switch dataType
        case 'Data'
            metadata.DataType = seal_struct.Type;
            metadata.DataSize = size(seal_struct.data);
            metadata.SamplingRate = seal_struct.srate;
            if isfield(seal_struct,'chanlocs')
                metadata.Chanlocs = 'Yes';
            end
        case 'Cortex'
            metadata.Cortex = 'Yes';
        case 'LeadField'
            metadata.LeadField = 'Yes';
            metadata.HeadModelType = seal_struct.HeadModelType;
            metadata.Orientation = seal_struct.Orientation;
        case 'Chanlocs'
            metadata.Chanlocs = 'Yes';
    end
        
    % --- Main Switch for different file formats ---
%     switch fileFormat
%         case 'eeglab_set'
%             if ~strcmp(dataType, 'Data')
%                 error('SEAL:ImportError', 'EEGLAB .set format is only applicable for EEG Data.');
%             end
%             seal_struct = import_from_eeglab(filePath);
%             
%         case 'bst_cortex'
%              if ~strcmp(dataType, 'Cortex')
%                 error('SEAL:ImportError', 'Brainstorm surface format is only applicable for Cortex data.');
%             end
%             seal_struct = import_from_bst_surface(filePath);
% 
%         case 'generic_mat'
%             sourceData = load(filePath); sourceData = sourceData{1};
%             seal_struct = map_from_generic_mat(sourceData, mappings);
%             
%         otherwise
%             error('SEAL:ImportError', 'Unsupported file format: %s', fileFormat);
%     end
end

% --- Private Helper Functions ---

function EEG_Data = import_from_eeglab(filePath)
    % Converts an EEGLAB .set file to the SEAL EEG_Data structure
    if ~exist('pop_loadset', 'file')
        error('SEAL:DependencyError', 'EEGLAB function pop_loadset not found. Is EEGLAB in your MATLAB path?');
    end
    EEG = pop_loadset('filename', filePath);
    
    EEG_Data = struct();
    EEG_Data.Comment = EEG.setname;
    EEG_Data.F = EEG.data;
    EEG_Data.srate = EEG.srate;
    EEG_Data.time = EEG.times / 1000; % Convert ms to seconds
    EEG_Data.chanlocs = EEG.chanlocs;
    EEG_Data.events = EEG.event;
    EEG_Data.type = 'EEG';
    EEG_Data.ica = struct();
end

function Cortex = import_from_bst_surface(filePath)
    % Converts a Brainstorm surface .mat file to the SEAL Cortex structure
    s = load(filePath);
    Cortex = struct();
    Cortex.Comment = s.Comment;
    Cortex.Vertices = s.Vertices;
    Cortex.Faces = s.Faces;
    if isfield(s, 'VertConn'), Cortex.VertConn = s.VertConn; end
end

function MappedStruct = map_from_generic_mat(sourceData, mappings, dataType)
    % Applies user-defined mappings to a generic struct
    MappedStruct = struct();
    allMapFields = fieldnames(mappings);

    % --- Process all provided mappings ---
    for i = 1:numel(allMapFields)
        sealField = allMapFields{i};
        mappingValue = mappings.(sealField);

        if ischar(mappingValue) && isfield(sourceData,mappingValue)
            % If the mapping is a string, treat it as a variable name (potentially nested)
            try
                MappedStruct.(sealField) = get_nested_field(sourceData, mappingValue);
            catch ME
                warning('SEAL:MappingWarning', 'Could not resolve source variable ''%s'' for field ''%s''. It will be skipped. Error: %s', mappingValue, sealField, ME.message);
            end
        else
            % If the mapping is not a string, assign the literal value
            MappedStruct.(sealField) = mappingValue;
        end
    end
    
    % --- Validate that all REQUIRED fields now exist in the MappedStruct ---
    requiredFields = {};
    switch dataType
        case 'Data', requiredFields = {'data', 'srate'};
        case 'Cortex', requiredFields = {'Vertices', 'Faces'};
        case 'LeadField', requiredFields = {'Gain', 'Orientation'};
    end
    
    for i = 1:numel(requiredFields)
        field = requiredFields{i};
        if ~isfield(MappedStruct, field)
            error('SEAL:MappingError', 'Required field ''%s'' was not successfully mapped and is missing from the final struct.', field);
        end
    end
    
    % --- Add Default Values ---
%     if ~isfield(MappedStruct, 'ID')
%         MappedStruct.ID = [dataType,'_default'];
%     end
end

function value = get_nested_field(s, fieldvalue)

    
    if isfield(s,fieldvalue)
        value = s.(fieldvalue);
    end
%     parts = strsplit(fieldPath, '.');
%     currentValue = s;
%     for k = 1:numel(parts)
%         if isfield(currentValue, parts{k})
%             currentValue = currentValue.(parts{k});
%         else
%             fprintf('Field path part ''%s'' not found.', parts{k});
%         end
%     end
%     value = currentValue;
>>>>>>> origin
=======
function [seal_struct, metadata] = seal_importer(filePath, dataType, varargin)
%SEAL_IMPORTER Core engine for importing data into SEAL format.
%   This function converts data from various formats into the standardized
%   SEAL Toolbox struct format. It is designed to be called directly from

%   scripts for batch processing.
%
%   Usage:
%   seal_struct = seal_importer_engine(filePath, fileFormat, dataType)
%   seal_struct = seal_importer_engine(filePath, 'generic_mat', dataType, mappings)
%
%   Inputs:
%   - filePath:   Full path to the source data file.
%   - fileFormat: String specifying the format (e.g., 'eeglab_set', 'bst_cortex', 'generic_mat').
%   - dataType:   The target SEAL data type ('Data', 'Cortex', 'LeadField').

%
%   Output:
%   - seal_struct: A standardized struct conforming to the SEAL data template.
metadata = struct('ID', [], 'DataType', [], 'DataSize', [], 'Chanlocs', [],...
                'SamplingRate', [], 'Cortex', [], 'LeadField', [], 'HeadModelType', [], 'Orientation', []);
Data_mapping = struct('ID', 'Data', 'data',[], 'Srate', [], 'Time', [], 'Chanlocs', [], 'Event', [], 'Type', 'EEG');
Cortex_mapping = struct('ID', 'Cortex  ', 'Vertices', [], 'Faces', [], 'VertConn', [], 'Atlas', [], 'Structure', [],'VertNormals',[], 'Type', 'Cortex');
LeadField_mapping = struct('ID', 'LeadField', 'Gain', [], 'Orientation', [],'GridLoc', [], 'GridOrient', [],...
    'HeadModelType', 'Surface', 'Type', 'LeadField');
Chanlocs_mapping = struct('ID', 'Chanlocs', 'chanlocs', [],'Type', 'Chanlocs');
default_mappings = struct('Data', Data_mapping, 'Cortex', Cortex_mapping, 'LeadField', LeadField_mapping,'Chanlocs',Chanlocs_mapping);

    % --- Input Validation ---
    p = inputParser;
    addRequired(p, 'filePath', @ischar);
%     addRequired(p, 'fileFormat', @ischar);
    addRequired(p, 'dataType', @ischar);
    addParameter(p, 'mappings', default_mappings.(dataType), @(x) isstruct(x));
    parse(p, filePath, dataType, varargin{:});
    
    mappings = p.Results.mappings;
    sourceData = load(filePath); 
    if numel(sourceData)==1
        sourceData = struct2cell(sourceData);
        sourceData = sourceData{1};
    end
    if ismatrix(sourceData)&&~isstruct(sourceData)
        sourceData = load(filePath);
    end
    seal_struct = map_from_generic_mat(sourceData, mappings , dataType);


    if ~isfield(seal_struct,'Type')
        seal_struct.Type = default_mappings.(dataType).Type;
    end
    switch dataType
        case 'Data'
            metadata.DataType = seal_struct.Type;
            metadata.DataSize = size(seal_struct.data);
            metadata.SamplingRate = seal_struct.srate;
            if isfield(seal_struct,'chanlocs')
                metadata.Chanlocs = 'Yes';
            end
        case 'Cortex'
            metadata.Cortex = 'Yes';
        case 'LeadField'
            metadata.LeadField = 'Yes';
            metadata.HeadModelType = seal_struct.HeadModelType;
            metadata.Orientation = seal_struct.Orientation;
        case 'Chanlocs'
            metadata.Chanlocs = 'Yes';
    end
        
    % --- Main Switch for different file formats ---
%     switch fileFormat
%         case 'eeglab_set'
%             if ~strcmp(dataType, 'Data')
%                 error('SEAL:ImportError', 'EEGLAB .set format is only applicable for EEG Data.');
%             end
%             seal_struct = import_from_eeglab(filePath);
%             
%         case 'bst_cortex'
%              if ~strcmp(dataType, 'Cortex')
%                 error('SEAL:ImportError', 'Brainstorm surface format is only applicable for Cortex data.');
%             end
%             seal_struct = import_from_bst_surface(filePath);
% 
%         case 'generic_mat'
%             sourceData = load(filePath); sourceData = sourceData{1};
%             seal_struct = map_from_generic_mat(sourceData, mappings);
%             
%         otherwise
%             error('SEAL:ImportError', 'Unsupported file format: %s', fileFormat);
%     end
end

% --- Private Helper Functions ---

function EEG_Data = import_from_eeglab(filePath)
    % Converts an EEGLAB .set file to the SEAL EEG_Data structure
    if ~exist('pop_loadset', 'file')
        error('SEAL:DependencyError', 'EEGLAB function pop_loadset not found. Is EEGLAB in your MATLAB path?');
    end
    EEG = pop_loadset('filename', filePath);
    
    EEG_Data = struct();
    EEG_Data.Comment = EEG.setname;
    EEG_Data.F = EEG.data;
    EEG_Data.srate = EEG.srate;
    EEG_Data.time = EEG.times / 1000; % Convert ms to seconds
    EEG_Data.chanlocs = EEG.chanlocs;
    EEG_Data.events = EEG.event;
    EEG_Data.type = 'EEG';
    EEG_Data.ica = struct();
end

function Cortex = import_from_bst_surface(filePath)
    % Converts a Brainstorm surface .mat file to the SEAL Cortex structure
    s = load(filePath);
    Cortex = struct();
    Cortex.Comment = s.Comment;
    Cortex.Vertices = s.Vertices;
    Cortex.Faces = s.Faces;
    if isfield(s, 'VertConn'), Cortex.VertConn = s.VertConn; end
end

function MappedStruct = map_from_generic_mat(sourceData, mappings, dataType)
    % Applies user-defined mappings to a generic struct
    MappedStruct = struct();
    allMapFields = fieldnames(mappings);

    % --- Process all provided mappings ---
    for i = 1:numel(allMapFields)
        sealField = allMapFields{i};
        mappingValue = mappings.(sealField);

        if ischar(mappingValue) && isfield(sourceData,mappingValue)
            % If the mapping is a string, treat it as a variable name (potentially nested)
            try
                MappedStruct.(sealField) = get_nested_field(sourceData, mappingValue);
            catch ME
                warning('SEAL:MappingWarning', 'Could not resolve source variable ''%s'' for field ''%s''. It will be skipped. Error: %s', mappingValue, sealField, ME.message);
            end
        else
            % If the mapping is not a string, assign the literal value
            MappedStruct.(sealField) = mappingValue;
        end
    end
    
    % --- Validate that all REQUIRED fields now exist in the MappedStruct ---
    requiredFields = {};
    switch dataType
        case 'Data', requiredFields = {'data', 'srate'};
        case 'Cortex', requiredFields = {'Vertices', 'Faces'};
        case 'LeadField', requiredFields = {'Gain', 'Orientation'};
    end
    
    for i = 1:numel(requiredFields)
        field = requiredFields{i};
        if ~isfield(MappedStruct, field)
            error('SEAL:MappingError', 'Required field ''%s'' was not successfully mapped and is missing from the final struct.', field);
        end
    end
    
    % --- Add Default Values ---
%     if ~isfield(MappedStruct, 'ID')
%         MappedStruct.ID = [dataType,'_default'];
%     end
end

function value = get_nested_field(s, fieldvalue)

    
    if isfield(s,fieldvalue)
        value = s.(fieldvalue);
    end
%     parts = strsplit(fieldPath, '.');
%     currentValue = s;
%     for k = 1:numel(parts)
%         if isfield(currentValue, parts{k})
%             currentValue = currentValue.(parts{k});
%         else
%             fprintf('Field path part ''%s'' not found.', parts{k});
%         end
%     end
%     value = currentValue;
>>>>>>> 3afa99beaec32453db712d4621fd05217f53fcfd
end
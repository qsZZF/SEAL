function defaultName = findAvailableName(basePath, namePattern)
%FINDAVAILABLENAME Find the smallest unused project name in a directory
%   DEFAULTNAME = FINDAVAILABLENAME(BASEPATH) searches in BASEPATH for
%   folders matching the pattern 'UNtitled N' and returns the smallest unused name.
%
%   defaultName = findAvailableName(BASEPATH, NAMEPATTERN) uses a custom
%   pattern where NAMEPATTERN should contain '#' as placeholder for numbers.
%
%   Examples:
%     findAvailableName('C:\Projects')
%     findAvailableName('/home/user/projects', 'MyProject_#')
%

    % Set default pattern if not provided
    if nargin < 2 || isempty(strtrim(namePattern))
        namePattern = 'Untitled #';
    end
    
    % Validate basePath
    if ~exist('basePath', 'var') || isempty(basePath) || ~exist(basePath, 'dir')
        defaultName = generateName(namePattern, 1);
        return;
    end
    
    % Extract the base name and prepare regex pattern
    [baseName, regexPattern] = preparePattern(namePattern);
    
    % Get all directories in the base path
    items = dir(basePath);
    dirNames = {items.name};
    
    % Remove '.' and '..'
    dirNames = dirNames(~ismember(dirNames, {'.', '..'}));
    
    if isempty(dirNames)
        defaultName = generateName(namePattern, 1);
        return;
    end
    
    % Find directories matching the pattern
    matchingDirs = {};
    numbers = [];
    
    for i = 1:length(dirNames)
        match = regexp(dirNames{i}, regexPattern, 'tokens', 'once');
        if ~isempty(match)
            matchingDirs{end+1} = dirNames{i};
            numbers(end+1) = round(str2double(match{2}));
        end
    end
    
    % Find the smallest available number
    if isempty(numbers)
        defaultName = generateName(namePattern, 1);
    else
        numbers = sort(numbers);
        nextNumber = findAvailableNumber(numbers);
        defaultName = generateName(namePattern, nextNumber);
    end
end

function [baseName, regexPattern] = preparePattern(namePattern)
%PREPAREPATTERN Convert name pattern to regex pattern
%   Replace '#' with regex digit pattern and extract base name
    
    if ~contains(namePattern, '#')
        % If no '#' found, use the pattern as base and append number
        baseName = namePattern;
        regexPattern = ['^', regexptranslate('escape', namePattern), ' (\d+)$'];
    else
        % Extract base name by replacing first '#' with empty
        baseName = strrep(namePattern, '#', '');
        
        % Create regex pattern by replacing '#' with digit capture group
        escapedPattern = regexptranslate('escape', namePattern);
        regexPattern = strrep(escapedPattern, '#', '(\d+)');
        regexPattern = ['^', regexPattern, '$'];
    end
end

function name = generateName(namePattern, number)
%GENERATENAME Generate project name from pattern and number
%   Replace '#' with number or append number if no '#' found
    if ~contains(namePattern, '#')
        name = sprintf('%s %d', namePattern, number);
    else
        name = strrep(namePattern, '#', num2str(number));
    end
end

function nextNumber = findAvailableNumber(numbers)
%FINDAVAILABLENUMBER Find the smallest available number in a sorted list
%   Handles gaps in numbering and finds the first available slot

    nextNumber = 1;
    for i = 1:length(numbers)
        if numbers(i) > nextNumber
            break;
        elseif numbers(i) == nextNumber
            nextNumber = nextNumber + 1;
        end
    end
end
function map = readPropertyFile(fileName)
%READPROPERTYFILE Read a set of properties from a file
%
%   PROPS = readPropertyFile(FILENAME)
%   Reads the set of properties stored in file FILENAME. Properties are
%   given in one line each, with a key, the equal sign, and the value of
%   the property. Comment in #-style are supported when they are aloe on
%   the line. The result is a Map object (implemented in containers.Map),
%   with each key being a property.
%   All properties are given as string.
%
%   Example
%     # This is a property file
%     prop1 = Some text
%     prop2 = 120
%     prop3 = [12 13 14]
%
%
%   Example
%   readPropertyFile
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-08-29,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

fID = fopen(fileName, 'wt');

map = containers.Map();

% Prints only the non-empty text lines
while true
    line = fgetl(fID);
    
    if line == -1
        break; 
    end
    
    line = strtrim(line);
    if isempty(line)
        continue; 
    end
    
    % check if line is a comment
    if line(1) == '#'
        continue;
    end
    
    % split each part of the line around the '=' sign
    tokens = regexp(line, '=', 'split');
    if length(tokens) ~= 2
        error('readPropertyFile:ParseError', ...
            ['Could not parse following line in file <' fileName '>:\n' line]);
    end
    
    key = strtrim(tokens{1});
    value = strtrim(tokens{2});
    
    if isKey(map, key)
        error('readPropertyFile:DuplicateKey', ...
            ['Property file <' fileName '> has duplicate key: ' key]);
    end
    
    map(key) = value;
end

fclose(fID);


function varargout = evaluate(obj, params)
%EVALUATE Evaluate the value and eventually the gradient of the function
%
%   output = evaluate(input)
%
%   Example
%   evaluate
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-01-06,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

% number of images
nImg = length(obj.Images);

% ensure we have initialized transformed images and gradients
if length(obj.TransformedImages)~=nImg
    createTransformedImages(obj);
end
if length(obj.TransformedGradients)~=nImg
    createTransformedGradients(obj);
end

% update transform parameters
ind1 = 1;
for i = 1:nImg
    % number of parameters of current transform
    transfo = obj.Transforms{i};
    nParams = getParameterLength(transfo);
    
    % extract parameters corresponding to current transform
    ind2 = ind1 + nParams-1;
    transfoParams = params(ind1:ind2);
    setParameters(transfo, transfoParams);
    
    % update for next transform
    ind1 = ind2 + 1;
end

% compute metric value and eventually gradient
if nargout <= 1
    varargout = {computeValue(obj)};
else
    [fval, grad] = computeValueAndGradient(obj);
    varargout = {fval, grad};
end

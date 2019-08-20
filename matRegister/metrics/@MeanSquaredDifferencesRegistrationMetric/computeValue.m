function [res, isInside] = computeValue(obj)
%COMPUTEVALUE Compute metric value.
%
% [VALUE, INSIDE] = METRIC.computeValue();
% Computes and return the value. Returns also a flag that indicates
% which test points belong to both images.
%

% compute values in image 1
[values1, inside1] = evaluate(obj.FixedImage, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.TransformedImage, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

% compute result
diff = (values2(isInside) - values1(isInside)).^2;
res = mean(diff);

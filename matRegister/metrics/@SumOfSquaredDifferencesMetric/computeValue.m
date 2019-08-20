function [res, isInside] = computeValue(obj)
%COMPUTEVALUE Compute metric value
%
% [VALUE, INSIDE] = computeValue(METRIC);
% Computes and return the value. Returns also a flag that indicates
% which test points belong to both images.
%

% compute values in image 1
[values1, inside1] = evaluate(obj.Img1, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.Img2, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

% compute result
% diff = (values1(isInside) - values2(isInside)).^2;
% res = sum(diff);
diff = (values1 - values2).^2;
res = sum(diff);

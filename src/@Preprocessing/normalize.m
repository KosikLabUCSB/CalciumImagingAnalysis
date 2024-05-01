function [P, yNorm] = normalize(obj, Y, varargin)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Preprocessing.normalize(X, varargin)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (October 2023)
% Maintainer: Ray Gifford up until May 2024

Y = Y - min(Y(:));      % remove negative offset

% data pre-processing
p = 1;
[P,yNorm] = preprocess_data(Y,p);

end
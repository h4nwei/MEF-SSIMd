function [Q,  qMap] = mef_ms_ssim_d(imgSeq, fI, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-exposure fused (MEF) for dynamic scene image quality model Version 1.0  %
% Copyright(c) 2019 Yuming Fang, Hanwei Zhu, Kede Ma, Zhou Wang, and Shutao Li  %
% All Rights Reserved.                                                          %
%                                                                               %
% ------------------------------------------------------------------------------%
% Permission to use, copy, or modify this software and its documentation        %
% for educational and research purposes only and without fee is hereby          %
% granted, provided that this copyright notice and the original authors'        %
% names appear on all copies and supporting documentation. This program         %
% shall not be used, rewritten, or adapted as the basis of a commercial         %
% software or hardware product without first obtaining permission of the        %
% authors. The authors make no representations about the suitability of         %
% this software for any purpose. It is provided "as is" without express         %
% or implied warranty.                                                          %
% ----------------------------------------------------------------------        %
% This is an implementation of an objective image quality assessment model      %
% for MEF of dynamic scene using their corresponding input source               %
% sequences as reference                                                        %
%                                                                               %
% Please refer to the following paper:                                          %
%                                                                               %
% Y. Fang, H. Zhu, K. Ma, Z. Wang, and S. Li, "Perceptual Evaluation for        %
% Multi-Exposure Image Fusion of Dynamic Scenes" submitted to IEEE              %
% Transactions on  Image Processing                                             %
%                                                                               %
%                                                                               %
% Kindly report any suggestions or corrections to h4nwei.zhu@gmail.com,         %
% fa0001ng@e.ntu.edu.sg, k29ma@uwaterloo.ca, or zhouwang@ieee.org               %
%                                                                               %
%                                                                               %
% ----------------------------------------------------------------------        %
% MEF-SSIMd                                                                      %
% input: (1) imgSeq: image sequences at multiple exposure levels [0-255].       %
%        (2) fI: the MEF image being evaluated in [0-255] grayscale.            %
% optional input:                                                               %
%        (3) C: constant in the SSIM index formula (see the above               %
%            reference). defualt value: K = (0.03*255)^2                        %
%        (4) p: the exponent parameter.  default value p = 4;                   %
%        (5) structureThres: the structure consistent threshold.                %
%                      defualt value: structureThres = 0.5                      %
%        (6) window: local window for statistics. default widnow is             %
%            Gaussian given by window = fspecial('gaussian', 11, 1.5);          %
%                                                                               %
%                                                                               %
%                                                                               %
% output: (1) oQ: The overlll quality score of the MEF image.                   %
%         (2)  Q: The quality scores in each scale.                             %
%         (3) qMap: The quality maps of the MEF image in each scale.            %
%                                                                               %
%Basic Usage:                                                                   %
%   Given the test MEF image and its corresponding source sequence              %
%                                                                               %
%   [Q, qMap] = mef_ms_ssim_d(imgSeq, fI);                                      %
%                                                                               %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input params parsing
params = inputParser;

default_C = (0.03*255)^2;
default_p = 4;
default_structureThres = 0.5;
default_window = fspecial('gaussian', 11, 1.5);

addRequired(params,'imgSeq');
addRequired(params,'fI');
addParameter(params, 'C', default_C, @isnumeric);
addParameter(params, 'p', default_p, @isnumeric);
addParameter(params, 'structureThres', default_structureThres, @isnumeric);
addParameter(params, 'window', default_window, @isnumeric);

parse(params, imgSeq, fI, varargin{:});

% initialization
C = params.Results.C;
p = params.Results.p;
structureThres = params.Results.structureThres;
window = params.Results.window;


imgSeq = double(imgSeq);
fI = double(fI);

% MEF-SSIMd
[Q, qMap] = mef_ssim_d(imgSeq, fI,  C, p, window, structureThres);

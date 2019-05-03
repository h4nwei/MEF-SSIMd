function [oQ,  qMap] = mef_ms_ssim_d(imgSeq, fI, varargin)
% ========================================================================

% input params parsing
params = inputParser;

default_C = (0.03*255)^2;
default_p = 4;
default_structureThres = 0.5;
default_weight = [0.0448  0.2856  0.3001]' / sum([0.0448  0.2856  0.3001]); 
default_window = fspecial('gaussian', 11, 1.5);
default_level = 1;
 
addRequired(params,'imgSeq');
addRequired(params,'fI');
addParameter(params, 'C', default_C, @isnumeric);
addParameter(params, 'p', default_p, @isnumeric);
addParameter(params, 'structureThres', default_structureThres, @isnumeric);
addParameter(params, 'weight', default_weight, @isnumeric);
addParameter(params, 'window', default_window, @isnumeric);
addParameter(params, 'level', default_level, @isnumeric);

parse(params, imgSeq, fI, varargin{:});

% initialization
C = params.Results.C;
p = params.Results.p;
structureThres = params.Results.structureThres;
weight = params.Results.weight;
window = params.Results.window;
level = params.Results.level;

[H, W] = size(window);

[s1, s2, s3] = size(imgSeq);
minImgWidth = min(s1, s2)/(2^(level-1));
maxWinWidth = max(H, W);
if (minImgWidth < maxWinWidth)
       return;
end

imgSeq = double(imgSeq);
fI = double(fI);
downsampleFilter = ones(2)./4;
Q1 = zeros(level,1);
Q2 = zeros(level,1);
qMap = cell(level,1);



if level == 1
    [Q1, qMap] = mef_ssim_d(imgSeq, fI,  C, p, window, structureThres);
    oQ = Q1;

    return;
else
    for l = 1 : level - 1
        [Q1(l),   qMap{l}] = mef_ssim_d(imgSeq, fI,  C, p, window, structureThres);
        imgSeqC = imgSeq;
        clear imgSeq;
        for i = 1:s3
            rI = squeeze(imgSeqC(:,:,i));
            dI = imfilter(rI, downsampleFilter, 'symmetric', 'same');
            imgSeq(:,:,i) = dI(1:2:end, 1:2:end);
        end
        dI = imfilter(fI, downsampleFilter, 'symmetric', 'same');
        clear fI;
        fI = dI(1:2:end, 1:2:end);
    end
    % the coarsest scale
    [Q1(level), qMap{level}] = mef_ssim_d(imgSeq, fI,  C, p, window, structureThres);
    Q1 = Q1(:);
    oQ = prod(Q1.^weight);  
end


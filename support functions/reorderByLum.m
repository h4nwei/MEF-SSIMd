function imgSeq_reordered = reorderByLum(imgSeq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function reorder each exposed image according to the expose time                          %
%   input:  1. imgSeq:  input image sequences at multiple exposure levels [0-255]               %
%                                                                                               %
%   output struct:                                                                              %
%           1. imgSeq_reorder: reordered image sequence at multiple exposure levels [0-255]     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgSeq = double(imgSeq);
m = squeeze(sum(sum(sum(imgSeq, 1), 2), 3));
[~, idx] = sort(m);
imgSeq_reordered = imgSeq(:, :, :, idx);
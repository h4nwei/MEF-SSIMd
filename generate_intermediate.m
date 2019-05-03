function [out_params] = generate_intermediate(imgSeqExd, C, p, window, structureThres, refIdx)

imgSeqExd = double(imgSeqExd);
[s1, s2, numExd] = size(imgSeqExd);
s3 = (numExd + 1) / 2;    % exposure number
wSize = size(window,1);
xIdxMax = s1-wSize+1;
yIdxMax = s2-wSize+1;
% C = (C*255)^2;
gMu = zeros(xIdxMax, yIdxMax, numExd); % global mean intensity
for i = 1 : numExd
    img = imgSeqExd(:,:,i);
    gMu(:,:,i) = ones(xIdxMax, yIdxMax) * mean(img(:)); %global mean intensity
end

lMu   = zeros(xIdxMax, yIdxMax, numExd); % local mean intensity
lMuSq = zeros(xIdxMax, yIdxMax, numExd);
sigmaSq = zeros(xIdxMax, yIdxMax, numExd); % signal strength from variance
for i = 1 : numExd
    lMu(:,:,i) = filter2(window, imgSeqExd(:, :, i), 'valid');
    lMuSq(:,:,i) = lMu(:,:,i) .* lMu(:,:,i); % mean square
    sigmaSq(:,:,i) = filter2(window, imgSeqExd(:, :, i).*imgSeqExd(:, :, i), 'valid') - lMuSq(:,:,i);
end

sigma = sqrt( max( sigmaSq, 0 ) );
ed = sigma*wSize  + 0.001; % signal strengh

% computing structural consistency map
count = 0;
numIndex = s3*(s3-1)/2;
sMap = zeros(xIdxMax, yIdxMax, s3, s3);
sRefMap = zeros(xIdxMax, yIdxMax, numIndex);
for i = 1 : s3
    for j = i+1 : s3
        count=count+1;
        crossMu = lMu(:,:,i) .* lMu(:,:,j);
        crossSigma = conv2(imgSeqExd(:, :, i).*imgSeqExd(:, :, j), window, 'valid') - crossMu;
        sMap(:,:,i,j) = ( crossSigma + C) ./ (sigma(:,:,i).* sigma(:,:,j) + C); % the third term in SSIM
        sRefMap(:,:,count) = sMap(:,:,i,j); % the third term in SSIM
    end
end
% sMap(sMap < 0) = 0;

% sRefMap = squeeze(sMap(:,:,refIdx,:)) + sMap(:,:,:,refIdx);
% sRefMap(:,:,refIdx) = ones(xIdxMax, yIdxMax); % add reference
sRefMap(sRefMap < structureThres) = 0;
sRefMap(sRefMap >= structureThres) = 1;
indexMap=zeros(xIdxMax,yIdxMax);

se = strel('disk',wSize);
for k = 1 : numIndex
    sRefMap(:,:,k) = imopen(sRefMap(:,:,k),se);
end

for k=1:numIndex
    indexMap=indexMap+sRefMap(:,:,k);
end
indexMap(indexMap <numIndex) = 0;
indexMap(indexMap >= numIndex) = 1;%final binary map
clear sMap;
clear sRefMap;
cMap = repmat(indexMap,[1, 1, s3]);
cMap(:,:,refIdx) = ones(xIdxMax, yIdxMax); % add reference
cMapExd = zeros(xIdxMax, yIdxMax, numExd);
cMapExd(:, :, 1:s3) = cMap;
clear cMap;
count = 0;
for i = 1 : s3
    if i ~= refIdx
        count = count + 1;
        cMapExd(:, :, count+s3) = 1 - cMapExd(:,:,i);
    end
end

% computing index matrix for the main loop
indM = zeros(xIdxMax, yIdxMax, s3);
indM(:,:,refIdx) = refIdx;
for i = 1 : s3
    if i < refIdx
        indM(:,:,i) = cMapExd(:, :, i) * i + cMapExd(:, :, i+s3) * (i+s3);
    elseif i > refIdx
        indM(:,:,i) = cMapExd(:, :, i) * i + cMapExd(:, :, i+s3-1) * (i+s3-1);
    end
end

% stepSize = 1;
% xIdx = 1 : stepSize : xIdxMax;
% xIdx = [xIdx xIdx(end)+1 : xIdxMax];
% yIdx = 1 : stepSize : yIdxMax;
% yIdx = [yIdx yIdx(end)+1 : yIdxMax];
% 
% offset = wSize-1;
% 
% R = zeros(xIdxMax, yIdxMax);
% for row = 1 : length(xIdx)
%     for col = 1 : length(yIdx)
%         i = xIdx(row);
%         j = yIdx(col);
%         blocks = imgSeqExd(i:i+offset, j:j+offset, indM(i,j,:));
%         vecs = reshape(blocks, [wSize*wSize, s3]);
%         denominator = 0;
%         for k=1:s3
%             denominator = denominator + norm(vecs(:,k) - lMu(i,j,k));
%         end
%         numerator = norm(sum(vecs,2) - mean(sum(vecs,2)));
%         R(i, j) = (numerator + eps) / (denominator + eps);
%     end
% end
% R(R > 1) = 1 - eps;
% R(R < 0) = 1 + eps;
% 
% w = tan(pi/2 * R);
% w(w > 10) = 10; %large number such as 10 is equivalent to taking maximum
% w = repmat(w, [1,1,numExd]);





sMap = ed.^p; % signal structure weighting map
sMap = sMap .* cMapExd + 0.001;
normalizer = sum(sMap,3);
sMap = sMap ./ repmat(normalizer,[1, 1, numExd]);

maxEd = ed .* cMapExd; %  desired signal strength
[maxEd, patchIndex] = max(maxEd, [], 3);

out_params.ed = ed;
out_params.lMu = lMu;
out_params.patchIndex = patchIndex;
out_params.sMap = sMap;
out_params.maxEd = maxEd;
out_params.indexMap = indexMap;
out_params.indM = indM;



